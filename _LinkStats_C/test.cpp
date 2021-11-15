#include "LinkStats.hpp"
#include <unistd.h>
#include <stdio.h>

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
    args.samFileName = cliArgs[1];
    args.fastaReferenceFileName = numCLIArgs > 2 ? cliArgs[2] : 0;
    args.overrideName = 0;
    args.fallbackName = "FallBackName";
    args.useMI = 1;
    args.arena = &workingSet;

    link_stats_return_data data;
    if (LinkStats(&args, data))
    {
        printf("\nGenome Length: %" PRIu64 "\n\n", data.genomeLength);
        for (const auto& [id, stats] : data.basicStats)
        {
            printf("%s:\n", id.c_str());
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
