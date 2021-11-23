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

#include "memory.hpp"

#include <stdlib.h>
#include <stdio.h>

global_function
u64
GetAlignmentPadding(u64 base, u32 alignment_pow2)
{
	u64 alignment = (u64)Pow2(alignment_pow2);
	u64 result = ((base + alignment - 1) & ~(alignment - 1)) - base;

	return(result);
}

void
CreateMemoryArena_(memory_arena *arena, u64 size, u32 alignment_pow2)
{
	u64 linkSize = sizeof(memory_arena);
	linkSize += GetAlignmentPadding(linkSize, alignment_pow2);
	u64 realSize = size + linkSize;

#ifndef _WIN32
	(void)posix_memalign((void **)&arena->base, Pow2(alignment_pow2), realSize);
#else
#include <memoryapi.h>
	(void)alignment_pow2;
	arena->base = (u08 *)VirtualAlloc(NULL, realSize, MEM_COMMIT | MEM_RESERVE, PAGE_READWRITE);
#endif
	arena->currentSize = 0;
	arena->maxSize = size;
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"	
	arena->next = (memory_arena *)arena->base;
#pragma clang diagnostic pop
	arena->base += linkSize;

	arena->next->base = 0;
	arena->active = 1;
}

void
FreeMemoryArena_(memory_arena *arena)
{
	if (arena->next)
	{
		if (arena->next->base)
		{
			FreeMemoryArena_(arena->next);
		}
		free(arena->next);
	}
}

void *
PushSize_(memory_arena *arena, u64 size, u32 alignment_pow2)
{
	if (!arena->active && arena->next && arena->next->base && !arena->next->currentSize)
	{
		arena->active = 1;
	}

	u64 padding = GetAlignmentPadding((u64)(arena->base + arena->currentSize), alignment_pow2);

	void *result;
	if (!arena->active || ((size + arena->currentSize + padding + sizeof(u64)) > arena->maxSize))
	{
		arena->active = 0;
		if (arena->next)
		{
			if (arena->next->base)
			{
				result = PushSize_(arena->next, size, alignment_pow2);
			}
			else
			{
				u64 linkSize = sizeof(memory_arena);
				linkSize += GetAlignmentPadding(linkSize, alignment_pow2);
				u64 realSize = size + padding + sizeof(u64) + linkSize;
				realSize = Max(realSize, arena->maxSize);

				CreateMemoryArenaP(arena->next, realSize, alignment_pow2);
				result = PushSize_(arena->next, size, alignment_pow2);
			}
		}
		else
		{
			result = 0;
			fprintf(stderr, "Push of %" PRIu64 " bytes failed, out of memory.\n", size);
			*((volatile u32 *)0) = 0;
		}
	}
	else
	{
		result = arena->base + arena->currentSize + padding;
		arena->currentSize += (size + padding + sizeof(u64));
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"		
		*((u64 *)(arena->base + arena->currentSize - sizeof(u64))) = (size + padding);
#pragma clang diagnostic pop
	}

	return(result);
}

void
FreeLastPush_(memory_arena *arena)
{
	if (!arena->active && arena->next && arena->next->base)
	{
		if (arena->next->active && !arena->next->currentSize)
		{
			arena->active = 1;
			FreeLastPush_(arena);
		}
		else
		{
			FreeLastPush_(arena->next);
		}
	}
	else if (arena->currentSize)
	{
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wcast-align"
		u64 sizeToRemove = *((u64 *)(arena->base + arena->currentSize - sizeof(u64)));
#pragma clang diagnostic pop		
		arena->currentSize -= (sizeToRemove + sizeof(u64));
	}
}

