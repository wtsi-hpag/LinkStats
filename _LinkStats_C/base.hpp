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

#include <inttypes.h>

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

#define global_function static
#define global_variable static

#define Min(x, y) ((x) < (y) ? (x) : (y))
#define Max(x, y) ((x) > (y) ? (x) : (y))

#define ArrayCount(array) (sizeof(array) / sizeof(array[0]))
#define ForLoop(n) for (u32 index = 0; index < (n); ++index)
#define TraverseLinkedList(startNode) u32 index = 0; for (auto *node = (startNode); node; node = node->next, ++index)
#define TraverseLinkedList2(startNode) u32 index2 = 0; for (auto *node2 = (startNode); node2; node2 = node2->next, ++index2)

#define Pow2(N) (1 << N)
