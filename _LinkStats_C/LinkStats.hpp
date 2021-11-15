#pragma once

#include <inttypes.h>

#include <map>
#include <string>
#include <vector>

typedef int8_t s08;
typedef int16_t s16;
typedef int32_t s32;
typedef int64_t s64;

typedef uint8_t u08;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef float f32;
typedef double f64;

#define KiloByte(x) (1024*x)
#define MegaByte(x) (1024*KiloByte(x))
#define GigaByte(x) (1024*MegaByte(x))

struct
memory_arena
{
   memory_arena *next;
   u08 *base;
   u64 currentSize;
   u64 maxSize;
   u64 active;
};

#define Default_Memory_Alignment_Pow2 4
void CreateMemoryArena_(memory_arena *arena, u64 size, u32 alignment_pow2 = Default_Memory_Alignment_Pow2);
#define CreateMemoryArena(arena, size, ...) CreateMemoryArena_(&arena, size, ##__VA_ARGS__)
#define CreateMemoryArenaP(arena, size, ...) CreateMemoryArena_(arena, size, ##__VA_ARGS__)

void FreeMemoryArena_(memory_arena *arena);
#define FreeMemoryArena(arena) FreeMemoryArena_(&arena)
#define FreeMemoryArenaP(arena) FreeMemoryArena_(arena)

struct
ll_node
{
    union
    {
        u64 value;
        void *ptr;
    };
    ll_node *next;
};

struct
ll
{
    ll_node *head;
    ll_node *tail;
    u64 count;
};

struct
basic_stats
{
    ll insertSizes;
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
    const char *samFileName;
    const char *fastaReferenceFileName;
    const char *overrideName;
    const char *fallbackName;
    memory_arena *arena;
    u08 useMI : 1;
};

struct
link_stats_return_data
{
    u64 genomeLength;
    std::vector<std::string> refNames;
    std::map<std::string, basic_stats *> basicStats;
    std::map<std::string, std::map<s32, std::map<std::pair<std::string, s32>, ll *>>> moleculeData;
};

u08 LinkStats(link_stats_run_args *args, link_stats_return_data &data);


