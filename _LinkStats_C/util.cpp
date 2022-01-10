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

#include "util.hpp"

#include <string.h>

bool
operator
<(u64_string const& lhs, u64_string const& rhs)
{
    if (lhs.id < rhs.id) return true;
    else if (lhs.id > rhs.id) return false;
    else
    {
        u64 *data1 = (u64 *)&lhs.string;
        u64 *data2 = (u64 *)&rhs.string;
        ForLoop(Min(lhs.length, rhs.length))
        {
            if (data1[index] < data2[index]) return true;
            else if (data1[index] > data2[index]) return false;
        }
        if (lhs.length < rhs.length) return true;
        else return false;
    }
}

bool
operator
==(u64_string const& lhs, u64_string const& rhs)
{
    if (lhs.id != rhs.id) return false;
    else if (lhs.length != rhs.length) return false;
    else
    {
        u64 *data1 = (u64 *)&lhs.string;
        u64 *data2 = (u64 *)&rhs.string;
        ForLoop(lhs.length) if (data1[index] != data2[index]) return false;
        return true;
    }
}

u64_string *
PushU64String(char *charString, memory_arena *arena, s32 id)
{
    if (!charString) return 0;

    u32 len = (((u32)strlen(charString)) >> 3) + 1;
    u64 *stringData = PushArrayP(arena, u64, len + 1);
    stringData[len] = 0;
    strcpy((char *)(stringData + 1), (const char *)charString);

    ((u32 *)stringData)[1] = len;
    ((s32 *)stringData)[0] = id;

    return (u64_string *)stringData;
}

char *
charU64String(u64_string *string)
{
    return (char *)&string->string;
}

void
MakeCopy(memory_arena *arena, u64_string *string, u64_string **out)
{
   *out = string;
}

bool
operator
<(char_string const& lhs, u64_string const& rhs)
{
    if (lhs.id < rhs.id) return true;
    else if (lhs.id > rhs.id) return false;
    else return strcmp((const char *)lhs.string, (const char *)charU64String((u64_string *)&rhs)) < 0;
}

bool
operator
==(char_string const& lhs, u64_string const& rhs)
{
    if (lhs.id != rhs.id) return false;
    else return !strcmp((const char *)lhs.string, (const char *)charU64String((u64_string *)&rhs));
}

void
MakeCopy(memory_arena *arena, char_string *string, u64_string **out)
{
    *out = PushU64String(string->string, arena, string->id);
}

void
MakeCopy(memory_arena *arena, s32 *id, s32 **out)
{
    s32 *tmp = PushArrayP(arena, s32, 1);
    tmp[0] = *id;
    *out = tmp;
}

insertsize_histogram *
NewInsertSizeHistogram(memory_arena *arena)
{
    insertsize_histogram *hist = PushStructP(arena, insertsize_histogram);
    memset(hist, 0, sizeof(*hist));
    
    return hist;
}

u64
EstimateMedian(insertsize_histogram *hist)
{
    s64 count = 0;
    ForLoop(InsertSize_Median_Estimate_Histogram_Size) count += (s64)hist->hist[index];
    count /= 2;
    u64 result;
    ForLoop(InsertSize_Median_Estimate_Histogram_Size)
    {
        result = (u64)index;
        count -= (s64)hist->hist[index];
        if (count <= 0) break;
    }

    return result;
}

bit_array *
CreateBitArray(u64 size, memory_arena *arena)
{
   bit_array *array = PushStructP(arena, bit_array);
   array->size = size;
   
   u64 nBytes = (size + 7) >> 3;
   array->bits = PushArrayP(arena, u08, nBytes);
   memset(array->bits, 0, nBytes);
   
   return array;
}

void
FillBitArray(bit_array *array, u64 from, u64 to)
{
   for (    u64 idx = from;
            idx < Min(to, array->size);
            ++idx )
      array->bits[idx >> 3] |= (1 << (idx & 7));
}

u08
IsBitSet(bit_array *array, u64 pos)
{
   return array->bits[pos >> 3] & (1 << (pos & 7));
}

ll<u64> *
FindGaps(bit_array *array, memory_arena *arena)
{
   ll<u64> *list;
   NewLL(arena, &list);

   u08 state1 = IsBitSet(array, 0);
   u64 len = 1;
   for (    u64 idx = 1;
            idx < array->size;
            ++idx )
   {
      u08 state2 = IsBitSet(array, idx);

      if (state1 && !state2) len = 0;
      else if (!state1 && state2) LLAddValue(list, len, arena);
      
      state1 = state2;
      ++len;
   }
   if (!state1) LLAddValue(list, len, arena);

   return list;
}
