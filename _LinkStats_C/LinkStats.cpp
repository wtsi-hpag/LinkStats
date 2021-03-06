/*
Copyright (c) 2022 Ed Harry, Wellcome Sanger Institute, Genome Research Limited

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

#include "LinkStats.hpp"

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wextra-semi-stmt"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconditional-uninitialized"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#define STB_SPRINTF_IMPLEMENTATION
#include "stb_sprintf.h"
#pragma clang diagnostic pop

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <htslib/kstring.h>

#include <stdlib.h>
#include <pthread.h>
#include <unistd.h>
#include <errno.h>

#define CircularBufferSize Pow2(10)
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

#define SAM_Data_Not_Sorted ((void *)1)
#define SAM_Data_Read_Error ((void *)2)

global_function
void *
FillCBTask(void *in)
{
    cb_fill_thread_data *data = (cb_fill_thread_data *)in;
    
    u32 bHead = 0;
    u32 head = (u32)data->buffer->head;
    
    u64 currPos = 0;
    s32 currTid = 0;
    u08 isSorted = 1;

    s32 readState;
    do
    {
        u32 headPlusOne = (head + 1) & (CircularBufferSize - 1);
        while (headPlusOne == (u32)data->buffer->tail) {}

        bam1_t *record = data->buffer->records + head;
        bam_set_mempolicy(record, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);

        if ((BufferSize - bHead) < AlignmentMemory)
        {
            while (!data->buffer->bufferTail) {}
            bHead = 0;
        }
        u32 freeMem;
        do
        {
            u32 bTail = (u32)data->buffer->bufferTail;
            freeMem = (((bHead < bTail) ? (bTail - 1) : BufferSize) - bHead) & (~7U);
        } while (freeMem < AlignmentMemory);

        record->data = data->buffer->dataBuffer + bHead;
        record->m_data = freeMem;

        if ((readState = sam_read1(data->samFile, data->header, record)) >= 0)
        {
            if (bam_get_mempolicy(record) & BAM_USER_OWNS_DATA) bHead += (record->m_data = (((u32)record->l_data + 7) & (~7U)));
            data->buffer->head = (volatile u32)(head = headPlusOne);

            s32 recTid;
            if ((recTid = (s32)record->core.tid) >= 0)
            {
                u64 recPos = (u64)record->core.pos;
                isSorted = (recTid > currTid) || ((recTid == currTid) && (recPos >= currPos)); 
                currPos = recPos;
                currTid = recTid;
            }
        }
    } while (isSorted && (readState >= 0));

    data->threadRunning = 0;
    pthread_exit(isSorted ? (readState < -1 ? SAM_Data_Read_Error : 0) : SAM_Data_Not_Sorted);
}

global_function
basic_stats *
NewBasicStats(memory_arena *arena)
{
    basic_stats *stats = PushStructP(arena, basic_stats);
    stats->medianInsertSize = 0;
    stats->totalReadLength = 0;
    stats->totalAlignments = 0;
    stats->totalDup = 0;
    stats->totalQCF = 0;
    stats->totalUnM = 0;
    stats->totalNoMI = 0;
    stats->totalNoBX = 0;
    stats->totalZeroMQ = 0;

    return stats;
}

global_function
void
UpdateMolecule(u64 groupCutOffDis, ll<molecule *> *molecules, bam1_t *record, s32 *mi, memory_arena *arena)
{
    molecule *mol = molecules->tail->data;
    u64 molMax = mol->maxCoord;
    
    u64 recMin = (u64)record->core.pos;
    u64 recMax = (u64)bam_endpos(record);
    u32 recReadLen = (u32)record->core.l_qseq; 
    u32 recQual = (u32)record->core.qual;

    u64 gap = (recMin > molMax) ? recMin - molMax : 0; 

    if (!molecules->count || (gap >= groupCutOffDis))
    {
        molecule *newMol = PushStructP(arena, molecule);
        newMol->minCoord = recMin;
        newMol->maxCoord = recMax;
        newMol->totalReadLength = recReadLen;
        newMol->nReads = 1;
        newMol->totalMappingQuality = recQual;
        if (mi)
        {
            newMol->haveMI = 1;
            newMol->mi = *mi;
        }
        else newMol->haveMI = 0;
        
        NewLL(arena, &newMol->gaps);

        LLAddValue(molecules, newMol, arena);
    }
    else
    {
        mol->minCoord = Min(mol->minCoord, recMin);
        mol->maxCoord = Max(molMax, recMax);
        mol->totalReadLength += recReadLen;
        mol->totalMappingQuality += recQual;
        if (mol->haveMI && (!mi || *mi != mol->mi)) mol->haveMI = 0;
        ++mol->nReads;

        if (gap) LLAddValue(mol->gaps, gap, arena);
    }
}

link_stats_return_data *
LinkStats(link_stats_run_args *args)
{
    if (!args) return 0;

    char LogBuffer[1024];
    u08 runState = 0;

    s32 logFd = args->logFD;
    u32 numThreads = Max(args->numThreads, 3);
    u64 groupCutOffDis = args->groupCutOffDis;
    u08 useMI = args->useMI;
    u64_string *overrideName = args->overrideName;
    u64_string *fallbackName = args->fallbackName;
    u64_string *fastaReferenceFileName = args->fastaReferenceFileName;
    u64_string *samFileName = args->samFileName;
    memory_arena *workingSet = args->arena;

    #define Log(message, ...) \
    do \
    { \
        u64 n = (u64)stbsp_snprintf(LogBuffer, sizeof(LogBuffer), message, ##__VA_ARGS__); \
        (void)write(logFd, LogBuffer, n); \
        (void)write(logFd, "\n", 1); \
    } while (0)

    if (!fallbackName)
    {
        Log("fallback name not set");
        return 0;
    }
    if (!samFileName)
    {
        Log("sam file name not set");
        return 0;
    }
    if (!workingSet)
    {
        Log("memory arena not set");
        return 0;
    }

    wavl_tree<u64_string, u64_string> *idLookup;
    InitialiseWavlTree(workingSet, &idLookup);

    wavl_tree<u64_string, insertsize_histogram> *insertsizeHists;
    InitialiseWavlTree(workingSet, &insertsizeHists);

    wavl_tree<u64_string, basic_stats> *basicStats; 
    InitialiseWavlTree(workingSet, &basicStats);

    wavl_tree<u64_string, wavl_tree<s32, wavl_tree<u64_string, ll<molecule *>>>> *moleculeData;
    InitialiseWavlTree(workingSet, &moleculeData);

    u64 genomeLength = 0;
    
    ll<u64_string *> *refNames;
    NewLL(workingSet, &refNames); 

    wavl_tree<u64_string, wavl_tree<s32, bit_array>> *coverage;
    InitialiseWavlTree(workingSet, &coverage);
    wavl_tree<u64_string, wavl_tree<s32, ll<u64>>> *coverageGaps;
    InitialiseWavlTree(workingSet, &coverageGaps);

    htsThreadPool threadPool = {0, (s32)(numThreads - 2)};
    htsFile *samFile = 0;
    sam_hdr_t *header = 0;
    if (    (threadPool.pool = hts_tpool_init(numThreads)) &&
            (samFile = hts_open(charU64String(samFileName), "r")) && 
            !hts_set_opt(samFile, HTS_OPT_THREAD_POOL, &threadPool) && 
            (!fastaReferenceFileName || !hts_set_opt(samFile, CRAM_OPT_REFERENCE, charU64String(fastaReferenceFileName))) &&
            (header = sam_hdr_read(samFile)))
    {
        kstring_t string = KS_INITIALIZE;
        if (!sam_hdr_find_tag_id(header, "HD", 0, 0, "SO", &string) && !strcmp((const char *)string.s, "coordinate"))
        {
            ForLoop((u32)header->n_targets) 
            {
                s32 idx = (s32)index;
                genomeLength += (u64)sam_hdr_tid2len(header, idx);
                LLAddValue(refNames, PushU64String((char *)sam_hdr_tid2name(header, idx), workingSet), workingSet);
            }

            circular_buffer *cBuffer = PushStructP(workingSet, circular_buffer);
            cBuffer->head = 0;
            cBuffer->tail = 0;
            cBuffer->bufferTail = 0;

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

                u64 recordCount = 0;
                u32 tail = (u32)cBuffer->tail;
                u32 head;
                while ((tail != (head = (u32)cBuffer->head)) || cbThreadData->threadRunning)
                {
                    if (tail != head)
                    {
                        bam1_t *record = cBuffer->records + tail;
                        
                        #define PassFilter(flags) (record->core.flag & (flags))
                        if (!PassFilter(BAM_FSUPPLEMENTARY | BAM_FSECONDARY))
                        {
                            u64_string *id;
                            if (!(id = overrideName))
                            {
                                char *tag;
                                u08 *tagData;
                                if ((tagData = bam_aux_get(record, "RG")) && (tag = bam_aux2Z(tagData)))
                                {
                                    char_string sTag(tag);
                                    auto *node = WavlTreeFindNode(workingSet, idLookup, &sTag);
                                    if (!node->value)
                                    {
                                        char *tmp;
                                        if (!sam_hdr_find_tag_id(header, "RG", "ID", tag, "SM", &string)) tmp = string.s;
                                        else tmp = tag;
                                        node->value = id = PushU64String(tmp, workingSet);
                                    }
                                    else id = node->value;
                                }
                                else id = fallbackName;
                            }

                            basic_stats *stats;
                            {
                                auto *node = WavlTreeFindNode(workingSet, basicStats, id);
                                if (!(stats = node->value)) node->value = stats = NewBasicStats(workingSet);
                            }
                           
                            insertsize_histogram *hist;
                            {
                                auto *node = WavlTreeFindNode(workingSet, insertsizeHists, id);
                                if (!(hist = node->value)) node->value = hist = NewInsertSizeHistogram(workingSet);
                            }
                            
                            s64 insertSize;
                            if (((insertSize = record->core.isize) > 0) && (insertSize < InsertSize_Median_Estimate_Histogram_Size)) ++hist->hist[insertSize];
                            
                            stats->totalReadLength += (u64)record->core.l_qseq;
                            ++stats->totalAlignments;

                            if (PassFilter(BAM_FUNMAP)) ++stats->totalUnM;
                            if (PassFilter(BAM_FDUP)) ++stats->totalDup;
                            if (PassFilter(BAM_FQCFAIL)) ++stats->totalQCF;

                            if (!PassFilter(BAM_FUNMAP | BAM_FDUP | BAM_FQCFAIL))
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
                                else
                                {
                                    bit_array *cov;
                                    {
                                        auto *node1 = WavlTreeFindNode(workingSet, coverage, id);
                                        if (!node1->value) InitialiseWavlTree(workingSet, &node1->value);
                                        auto *node2 = WavlTreeFindNode(workingSet, node1->value, &record->core.tid);
                                        if (!(cov = node2->value)) node2->value = cov = CreateBitArray((u64)sam_hdr_tid2len(header, record->core.tid), workingSet);
                                    }
                                    
                                    FillBitArray(cov, (u64)record->core.pos, (u64)bam_endpos(record));
                                }
                                
                                if (record->core.qual && bx && (useMI ? mi : (s32 *)1))
                                {
                                    auto *node1 = WavlTreeFindNode(workingSet, moleculeData, id);
                                    if (!node1->value) InitialiseWavlTree(workingSet, &node1->value); 
                                    auto *node2 = WavlTreeFindNode(workingSet, node1->value, &record->core.tid);
                                    if (!node2->value) InitialiseWavlTree(workingSet, &node2->value);
                                    char_string bxmi(bx, useMI ? *mi : 0);
                                    auto *node3 = WavlTreeFindNode(workingSet, node2->value, &bxmi);
                                    if (!node3->value) NewLL(workingSet, &node3->value);

                                    UpdateMolecule(groupCutOffDis, node3->value, record, mi, workingSet);
                                }
                            }
                        }

                        if (bam_get_mempolicy(record) & BAM_USER_OWNS_DATA) cBuffer->bufferTail = (volatile u32)((u32)(record->data - cBuffer->dataBuffer) + record->m_data);
                        else bam_destroy1(record);

                        cBuffer->tail = (volatile u32)(tail = ((tail + 1) & (CircularBufferSize - 1)));

                        {
                            #define Log2_Print_Interval 14
                            if (!(++recordCount & ((Pow2(Log2_Print_Interval)) - 1)))
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

                void *fileReadError = 0;
                if (pthread_join(readingThread, &fileReadError) || fileReadError)
                {
                    Log("Read error");
                    if (fileReadError == SAM_Data_Not_Sorted) goto NotSorted;
                }
                else
                {
                    runState = 1;

                    {
                        TraverseLinkedList(WavlTreeFreezeToLL(insertsizeHists))  WavlTreeFindNode(workingSet, basicStats, node->key)->value->medianInsertSize = EstimateMedian(node->value);
                    }
                    
                    {
                        TraverseLinkedList(WavlTreeFreezeToLL(coverage))
                        {
                            auto *node3 = WavlTreeFindNode(workingSet, coverageGaps, node->key);
                            InitialiseWavlTree(workingSet, &node3->value);
                            TraverseLinkedList2(WavlTreeFreezeToLL(node->value)) WavlTreeFindNode(workingSet, node3->value, node2->key)->value = FindGaps(node2->value, workingSet);
                        }
                    }

                    Log("\nGenome Length: %$" PRIu64 "bp\n", genomeLength);
                    TraverseLinkedList(WavlTreeFreezeToLL(basicStats))
                    {
                        basic_stats *stats = node->value;
                        Log("%s:", charU64String(node->key));
                        Log("\t%$" PRIu64 "bp median insert size", stats->medianInsertSize);
                        Log("\t%$" PRIu64 "bp total read length", stats->totalReadLength);
                        Log("\t%$" PRIu64 " alignments", stats->totalAlignments);
                        Log("\t%$" PRIu64 " unmapped", stats->totalUnM);
                        Log("\t%$" PRIu64 " duplicates", stats->totalDup);
                        Log("\t%$" PRIu64 " QC fail", stats->totalQCF);
                        Log("\t%$" PRIu64 " no BX", stats->totalNoBX);
                        Log("\t%$" PRIu64 " no MI", stats->totalNoMI);
                        Log("\t%$" PRIu64 " 0 mapq", stats->totalZeroMQ);
                        Log("");
                    }
                }
            }
        }
        else
        {
NotSorted:
            Log("Cannot process non-coordinate-sorted SAM data");
        }
        ks_free(&string);
    }

    if (header) sam_hdr_destroy(header);
    if (samFile && hts_close(samFile)) runState = 0;
    hts_tpool_destroy(threadPool.pool);
    
    link_stats_return_data *returnData;
    if ((returnData = (runState ? PushStructP(workingSet, link_stats_return_data) : 0)))
    {
        returnData->genomeLength = genomeLength;
        returnData->refNames = refNames;
        returnData->basicStats = basicStats;
        returnData->moleculeData = moleculeData;
        returnData->coverageGaps = coverageGaps;
    }

    return returnData;
}
