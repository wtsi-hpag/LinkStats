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

#pragma once

#include "base.hpp"

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

void *PushSize_(memory_arena *arena, u64 size, u32 alignment_pow2 = Default_Memory_Alignment_Pow2);
#define PushStruct(arena, type, ...) (type *)PushSize_(&arena, sizeof(type), ##__VA_ARGS__)
#define PushArray(arena, type, n, ...) (type *)PushSize_(&arena, sizeof(type) * n, ##__VA_ARGS__)
#define PushStructP(arena, type, ...) (type *)PushSize_(arena, sizeof(type), ##__VA_ARGS__)
#define PushArrayP(arena, type, n, ...) (type *)PushSize_(arena, sizeof(type) * n, ##__VA_ARGS__)

void FreeLastPush_(memory_arena *arena);
#define FreeLastPush(arena) FreeLastPush_(&arena)
#define FreeLastPushP(arena) FreeLastPush_(arena)
