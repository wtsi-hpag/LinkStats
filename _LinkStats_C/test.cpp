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

#include "LinkStats.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

s32
main(s32 numCLIArgs, const char **cliArgs)
{
    s32 exitCode = EXIT_FAILURE;
    
    if (numCLIArgs <= 1) return exitCode;

    memory_arena workingSet;
    CreateMemoryArena(workingSet, MegaByte(512));
    
    link_stats_run_args args;
    args.logFD = STDERR_FILENO;
    args.numThreads = 8;
    args.samFileName = PushU64String((char *)cliArgs[1], &workingSet);
    args.fastaReferenceFileName = numCLIArgs > 2 ? PushU64String((char *)cliArgs[2], &workingSet) : 0;
    args.overrideName = 0;
    args.fallbackName = PushU64String((char *)"FallBackName", &workingSet);
    args.useMI = 1;
    args.arena = &workingSet;

    link_stats_return_data data;
    if (LinkStats(&args, &data))
    {
        printf("\nGenome Length: %" PRIu64 "\n\n", data.genomeLength);
        TraverseLinkedList(WavlTreeFreezeToLL(data.basicStats))
        {
            basic_stats *stats = node->value;
            printf("%s:\n", charU64String(node->key));
            printf("\t%" PRIu64 " inserts\n", stats->insertSizes.count);
            printf("\t%" PRIu64 " total read length\n", stats->totalReadLength);
            printf("\t%" PRIu64 " alignments\n", stats->totalAlignments);
            printf("\t%" PRIu64 " unmapped\n", stats->totalUnM);
            printf("\t%" PRIu64 " duplicates\n", stats->totalDup);
            printf("\t%" PRIu64 " QC fail\n", stats->totalQCF);
            printf("\t%" PRIu64 " no BX\n", stats->totalNoBX);
            printf("\t%" PRIu64 " no MI\n", stats->totalNoMI);
            printf("\t%" PRIu64 " 0 mapq\n", stats->totalZeroMQ);
            printf("\n");    
        }

        exitCode = EXIT_SUCCESS;
    }

    FreeMemoryArena(workingSet);
    return exitCode;
}
