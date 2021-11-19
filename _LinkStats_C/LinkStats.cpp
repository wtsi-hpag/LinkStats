/*
Copyright (c) 2021 Ed Harry, Wellcome Sanger Institute, Genome Research Limited

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "Header.h"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/kstring.h>

#define CircularBufferSize (1 << 10)
#define AlignmentMemory 512
#define BufferSize (CircularBufferSize * AlignmentMemory)

struct
circular_buffer
{
    bam1_t records[CircularBufferSize];
    u08 dataBuffer[BufferSize];
    volatile u32 head;
    volatile u32 tail;
    volatile u32 bufferTail;
};

struct
cb_fill_thread_data
{
    htsFile *samFile;
    sam_hdr_t *header;
    circular_buffer *buffer;
    volatile u08 threadRunning;
};

global_function
void *
FillCBTask(void *in)
{
    cb_fill_thread_data *data = (cb_fill_thread_data *)in;
    u32 bHead = 0;
    u32 head = (u32)data->buffer->head;

    s32 readState;
    do
    {
	u32 headPlusOne = (head + 1) & (CircularBufferSize - 1);
	while (headPlusOne == (u32)data->buffer->tail) {}

	bam1_t *record = data->buffer->records + head;
	bam_set_mempolicy(record, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);
	
	if ((BufferSize - bHead) < AlignmentMemory) bHead = 0;
	u32 freeMem;
	do
	{
	    u32 bTail = (u32)data->buffer->bufferTail;
	    freeMem = (((bHead < bTail) ? bTail : BufferSize) - bHead) & (~7U);
	} while (freeMem < AlignmentMemory);
	
	record->data = data->buffer->dataBuffer + bHead;
	record->m_data = freeMem;

	if ((readState = sam_read1(data->samFile, data->header, record)) >= 0)
	{
	    if (bam_get_mempolicy(record) & BAM_USER_OWNS_DATA) bHead += (record->m_data = (((u32)record->l_data + 7) & (~7U)));
	    data->buffer->head = (volatile u32)(head = headPlusOne);
	}
    } while (readState >= 0);

    data->threadRunning = 0;
    pthread_exit((void *)(readState < -1));
}

global_function
u16
PassFilter(bam1_t *record, u16 flags)
{
    return record->core.flag & flags;
}

#define DeclareLLAddFunction(type, id) \
    global_function \
    void \
    LLAddValue(ll *list, type data, memory_arena *arena) \
{ \
    ll_node *newNode = PushStructP(arena, ll_node); \
    newNode->id = data; \
    newNode->next = 0; \
    \
    if (!list->head) list->head = newNode; \
    else list->tail->next = newNode; \
    list->tail = newNode; \
    \
    ++list->count; \
}\

DeclareLLAddFunction(u64, value)
DeclareLLAddFunction(void *, ptr)

global_function
basic_stats *
NewBasicStats(memory_arena *arena)
{
    basic_stats *stats = PushStructP(arena, basic_stats);
    memset(stats, 0, sizeof(basic_stats));
    return stats;
}

global_function
ll *
NewLL(memory_arena *arena)
{
    ll *list = PushStructP(arena, ll);
    memset(list, 0, sizeof(ll));
    return list;
}

global_function
alignment *
NewAlignment(bam1_t *record, s32 *mi, memory_arena *arena)
{
    alignment *al = PushStructP(arena, alignment);
    
    al->haveMI = mi ? 1 : 0;
    al->qual = (u32)record->core.qual;
    al->referenceStart = (u64)record->core.pos;
    al->referenceEnd = (u64)bam_endpos(record);
    al->queryLength = (u64)record->core.l_qseq;
    al->mi = mi ? *mi : 0;

    return al;
}

u08
LinkStats(link_stats_run_args *args, link_stats_return_data &returnData)
{
    char LogBuffer[1024];
    u08 runState = 0;

    s32 logFd = args->logFD;
    u32 numThreads = Max(args->numThreads, 3);
    u08 useMI = args->useMI;
    const char *overrideName = args->overrideName;
    const char *fallbackName = args->fallbackName;
    const char *fastaReferenceFileName = args->fastaReferenceFileName;
    const char *samFileName = args->samFileName;
    memory_arena *workingSet = args->arena;

    #define Log(message, ...) \
    do \
    { \
	u64 n = (u64)stbsp_snprintf(LogBuffer, sizeof(LogBuffer), message, ##__VA_ARGS__); \
	write(logFd, LogBuffer, n); \
	write(logFd, "\n", 1); \
    } while (0)

    std::map<std::string, std::string> idLookup;
    std::map<std::string, basic_stats *> basicStats;
    std::map<std::string, std::map<s32, std::map<std::pair<std::string, s32>, ll *>>> moleculeData;
    u64 genomeLength = 0;
    std::vector<std::string> refNames;

    htsThreadPool threadPool = {0, (s32)(numThreads - 2)};
    htsFile *samFile;
    sam_hdr_t *header;
    if (    (threadPool.pool = hts_tpool_init(numThreads)) &&
	    (samFile = hts_open(samFileName, "r")) && 
	    !hts_set_opt(samFile, HTS_OPT_THREAD_POOL, &threadPool) && 
	    (!fastaReferenceFileName || !hts_set_opt(samFile, CRAM_OPT_REFERENCE, fastaReferenceFileName)) &&
	    (header = sam_hdr_read(samFile)))
    {
	ForLoop((u32)header->n_targets) 
	{
	    genomeLength += (u64)sam_hdr_tid2len(header, (s32)index);
	    refNames.push_back(sam_hdr_tid2name(header, (s32)index));
	}
	
	circular_buffer *cBuffer = PushStructP(workingSet, circular_buffer);
	cBuffer->head = cBuffer->tail = cBuffer->bufferTail = 0;

	cb_fill_thread_data *cbThreadData = PushStructP(workingSet, cb_fill_thread_data);
	cbThreadData->samFile = samFile;
	cbThreadData->header = header;
	cbThreadData->buffer = cBuffer;
	cbThreadData->threadRunning = 1;

	pthread_t readingThread;
	if (!pthread_create(&readingThread, 0, FillCBTask, cbThreadData))
	{
	    char printNBuffers[2][16] = {{0}};
	    u08 printNBufferPtr = 0;

	    kstring_t string = KS_INITIALIZE;
	    u64 recordCount = 0;
	    u32 tail = (u32)cBuffer->tail;
	    u32 head;
	    while ((tail != (head = (u32)cBuffer->head)) || cbThreadData->threadRunning)
	    {
		if (tail != head)
		{
		    bam1_t *record = cBuffer->records + tail;

		    if (!PassFilter(record, BAM_FSUPPLEMENTARY | BAM_FSECONDARY))
		    {
			char *id;
			if (!(id = (char *)overrideName))
			{
			    char *tag;
			    u08 *tagData;
			    if ((tagData = bam_aux_get(record, "RG")) && (tag = bam_aux2Z(tagData)))
			    {
				if (idLookup.find(tag) == idLookup.end())
				{
				    if (!sam_hdr_find_tag_id(header, "RG", "ID", tag, "SM", &string)) id = string.s;
				    else id = tag;
				    idLookup[tag] = id;
				}
				else id = (char *)idLookup[tag].c_str();
			    }
			    else id = (char *)fallbackName;
			}

			basic_stats *stats;
			if (!(stats = basicStats[id])) basicStats[id] = stats = NewBasicStats(workingSet);

			{
			    s64 insertSize;
			    if ((insertSize = record->core.isize) > 0) LLAddValue(&stats->insertSizes, (u64)insertSize, workingSet);

			    stats->totalReadLength += (u64)record->core.l_qseq;
			    ++stats->totalAlignments;

			    if (PassFilter(record, BAM_FUNMAP)) ++stats->totalUnM;
			    if (PassFilter(record, BAM_FDUP)) ++stats->totalDup;
			    if (PassFilter(record, BAM_FQCFAIL)) ++stats->totalQCF;
			}

			if (!PassFilter(record, BAM_FUNMAP | BAM_FDUP | BAM_FQCFAIL))
			{
			    char *bx = 0;
			    {
				u08 *bxData;
				if ((bxData = bam_aux_get(record, "BX"))) bx = bam_aux2Z(bxData);
			    }
			    s32 *mi = 0;
			    s32 miValue;
			    {
				u08 *miData;
				if ((miData = bam_aux_get(record, "MI")))
				{
				    miValue = bam_aux2i(miData);
				    if (errno != EINVAL) mi = &miValue;
				}
			    }

			    if (!bx) ++stats->totalNoBX;
			    if (!mi) ++stats->totalNoMI;
			    if (!record->core.qual) ++stats->totalZeroMQ;

			    if (record->core.qual && bx && (useMI ? mi : (s32 *)1))
			    {
				alignment *al = NewAlignment(record, mi, workingSet);
				std::pair<std::string, s32> bxmi = std::make_pair(bx, useMI ? *mi : 0);

				ll *list;
				if (!(list = moleculeData[id][record->core.tid][bxmi])) moleculeData[id][record->core.tid][bxmi] = list = NewLL(workingSet);

				LLAddValue(list, al, workingSet);
			    }
			}
		    }

		    if (bam_get_mempolicy(record) & BAM_USER_OWNS_DATA) cBuffer->bufferTail = (volatile u32)((u32)(record->data - cBuffer->dataBuffer) + record->m_data);
		    else bam_destroy1(record);

		    cBuffer->tail = (volatile u32)(tail = ((tail + 1) & (CircularBufferSize - 1)));

		    {
			#define Log2_Print_Interval 14
			if (!(++recordCount & ((1 << Log2_Print_Interval) - 1)))
			{
			    u08 currPtr = printNBufferPtr;
			    u08 otherPtr = (currPtr + 1) & 1;
			    stbsp_snprintf(printNBuffers[currPtr], sizeof(printNBuffers[currPtr]), "%$" PRIu64, recordCount);

			    if (strcmp(printNBuffers[currPtr], printNBuffers[otherPtr])) Log("%s reads processed", printNBuffers[currPtr]);

			    printNBufferPtr = otherPtr;
			}
		    }
		}
	    }

	    void *fileReadError;
	    if (pthread_join(readingThread, &fileReadError) || fileReadError) Log("Read error");
	    else runState = 1;

	    {
		Log("\nGenome Length: %" PRIu64 "\n", genomeLength);
		for (const auto& [id, stats] : basicStats)
		{
		    Log("%s:", id.c_str());
		    Log("\t%" PRIu64 " inserts", stats->insertSizes.count);
		    Log("\t%" PRIu64 " total read length", stats->totalReadLength);
		    Log("\t%" PRIu64 " alignments", stats->totalAlignments);
		    Log("\t%" PRIu64 " unmapped", stats->totalUnM);
		    Log("\t%" PRIu64 " duplicates", stats->totalDup);
		    Log("\t%" PRIu64 " QC fail", stats->totalQCF);
		    Log("\t%" PRIu64 " no BX", stats->totalNoBX);
		    Log("\t%" PRIu64 " no MI", stats->totalNoMI);
		    Log("\t%" PRIu64 " 0 mapq", stats->totalZeroMQ);
		    Log("");
		}
	    }

	    ks_free(&string);
	}

	sam_hdr_destroy(header);
	if (hts_close(samFile)) runState = 0;
	hts_tpool_destroy(threadPool.pool);
    }

    if (runState)
    {
	returnData.genomeLength = genomeLength;
	returnData.refNames = refNames;
	returnData.basicStats = basicStats;
	returnData.moleculeData = moleculeData;
    }

    return runState;
}
