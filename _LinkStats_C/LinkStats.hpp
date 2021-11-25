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

#pragma once

#include "base.hpp"
#include "memory.hpp"
#include "util.hpp"

struct
basic_stats
{
    ll<u64> insertSizes;
    u64 totalReadLength;
    u64 totalAlignments;
    u64 totalDup;
    u64 totalQCF;
    u64 totalUnM;
    u64 totalNoMI;
    u64 totalNoBX;
    u64 totalZeroMQ;
};

struct
alignment
{
    u32 haveMI : 1;
    u32 qual : 31;
    u64 referenceStart;
    u64 referenceEnd;
    u64 queryLength;
    s32 mi;
};

struct
link_stats_run_args
{
    s32 logFD;
    u32 numThreads;
    u64_string *samFileName;
    u64_string *fastaReferenceFileName;
    u64_string *overrideName;
    u64_string *fallbackName;
    memory_arena *arena;
    u08 useMI : 1;
};

struct
link_stats_return_data
{
    u64 genomeLength;
    ll<u64_string *> *refNames;
    wavl_tree<u64_string, basic_stats> *basicStats;
    wavl_tree<u64_string, wavl_tree<s32, wavl_tree<u64_string, ll<alignment *>>>> *moleculeData;
};

link_stats_return_data *LinkStats(link_stats_run_args *args);

