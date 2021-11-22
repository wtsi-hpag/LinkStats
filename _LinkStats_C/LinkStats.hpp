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

#include <stdlib.h>
#include <inttypes.h>
#include <string.h>

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

void *PushSize_(memory_arena *arena, u64 size, u32 alignment_pow2 = Default_Memory_Alignment_Pow2);
#define PushStruct(arena, type, ...) (type *)PushSize_(&arena, sizeof(type), ##__VA_ARGS__)
#define PushArray(arena, type, n, ...) (type *)PushSize_(&arena, sizeof(type) * n, ##__VA_ARGS__)
#define PushStructP(arena, type, ...) (type *)PushSize_(arena, sizeof(type), ##__VA_ARGS__)
#define PushArrayP(arena, type, n, ...) (type *)PushSize_(arena, sizeof(type) * n, ##__VA_ARGS__)

struct
ll_node
{
    ll_node *next;
    union
    {
        u64 value;
        void *ptr;
    };
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
u64_string
{
    s32 id;
    u32 length;
    u64 string;
    friend bool operator< (u64_string const& a, u64_string const& b);
    friend bool operator== (u64_string const& a, u64_string const& b);
};


u64_string *PushU64String(char *charString, memory_arena *arena, s32 id = 0);
char *charU64String(u64_string *string);
void MakeCopy(memory_arena *arena, u64_string *string, u64_string **out);

struct
char_string
{
   char *string;
   s32 id;
   char_string(char *s, s32 i = 0) { string = s; id = i; }
   friend bool operator< (char_string const& a, char_string const& b);
   friend bool operator== (char_string const& a, char_string const& b);
};

void MakeCopy(memory_arena *arena, char_string *string, char_string **out);
void MakeCopy(memory_arena *arena, char_string *string, u64_string **out);
void MakeCopy(memory_arena *arena, s32 *tid, s32 **out);

template
<typename k, typename v>
struct
wavl_node
{
    wavl_node<k,v> *parent;
    wavl_node<k,v> *left;
    wavl_node<k,v> *right;
    u64 rank;
    k *key;
    v *value;
};

template
<typename k, typename v>
struct
wavl_tree
{
    wavl_node<k,v> *root;
    u64 size;
};

template
<typename k, typename v>
void
InitialiseWavlTree(memory_arena *arena, wavl_tree<k,v> **out)
{
    wavl_tree<k,v> *tree = (wavl_tree<k,v> *)PushSize_(arena, sizeof(wavl_tree<k,v>));
    tree->size = 0;
    tree->root = 0;
    
    *out = tree;
}

template
<typename k, typename v, typename l>
wavl_node<k,v> *
WavlTreeNewNode(memory_arena *arena, wavl_node<k,v> *parent, l *key)
{
    wavl_node<k,v> *node = (wavl_node<k,v> *)PushSize_(arena, sizeof(wavl_node<k,v>));
    node->parent = parent;
    node->left = 0;
    node->right = 0;
    node->value = 0;
    MakeCopy(arena, key, &node->key);
    node->rank = 0;

    return(node);
}

template
<typename k, typename v>
u32
WavlTreeNeedToRotateLeftStrong(wavl_node<k,v> *node)
{
    return((node->rank == 1 && !node->left) || (node->rank > node->left->rank + 1));
}

template
<typename k, typename v>
u32
WavlTreeNeedToRotateRightStrong(wavl_node<k,v> *node)
{
    return((node->rank == 1 && !node->right) || (node->rank > node->right->rank + 1));
}

template
<typename k, typename v>
u32
WavlTreeNeedToRotateLeftWeak(wavl_node<k,v> *node)
{
    return((!node->left) || (node->rank > node->left->rank + 1));
}

template
<typename k, typename v>
u32
WavlTreeNeedToRotateRightWeak(wavl_node<k,v> *node)
{
    return((!node->right) || (node->rank > node->right->rank + 1));
}

template
<typename k, typename v>
void
WavlTreeRotateLeft(wavl_tree<k,v> *tree, wavl_node<k,v> *node)
{
    wavl_node<k,v> *right = node->right;
    node->right = right->left;
    if (right->left)
    {
        right->left->parent = node;
    }
    right->parent = node->parent;
    if (!node->parent)
    {
        tree->root = right;
    }
    else if (node->parent->left == node)
    {
        node->parent->left = right;
    }
    else
    {
        node->parent->right = right;
    }
    right->left = node;
    node->parent = right;
}

template
<typename k, typename v>
void
WavlTreeRotateRight(wavl_tree<k,v> *tree, wavl_node<k,v> *node)
{
    wavl_node<k,v> *left = node->left;
    node->left = left->right;
    if (left->right)
    {
        left->right->parent = node;
    }
    left->parent = node->parent;
    if (!node->parent)
    {
        tree->root = left;
    }
    else if (node->parent->right == node)
    {
        node->parent->right = left;
    }
    else
    {
        node->parent->left = left;
    }
    left->right = node;
    node->parent = left;
}

template
<typename k, typename v>
void
WavlTreeBalance(wavl_tree<k,v> *tree, wavl_node<k,v> *node)
{
    for (   wavl_node<k,v> *parent = node->parent;
            parent && ((node->rank + 1) != parent->rank);
            node = parent, parent = node->parent, ++node->rank )
    {
        if (parent->left == node)
        {
            if (WavlTreeNeedToRotateRightStrong(parent))
            {
                if (WavlTreeNeedToRotateLeftWeak(node))
                {
                    --node->rank;
                    ++node->right->rank;
                    WavlTreeRotateLeft(tree, node);
                }
                --parent->rank;
                WavlTreeRotateRight(tree, parent);
                break;
            }
        }
        else if (WavlTreeNeedToRotateLeftStrong(parent))
        {
            if (WavlTreeNeedToRotateRightWeak(node))
            {
                --node->rank;
                ++node->left->rank;
                WavlTreeRotateRight(tree, node);
            }
            --parent->rank;
            WavlTreeRotateLeft(tree, parent);
            break;
        }
    }
}

template
<typename k, typename v, typename l>
wavl_node<k,v> *
WavlTreeInsertValue(memory_arena *arena, wavl_tree<k,v> *tree, l *key)
{
    wavl_node<k,v> *node = tree->root;
    if (!node)
    {
        ++tree->size;
        return (tree->root = WavlTreeNewNode(arena, (wavl_node<k,v> *)0, key));
    }
    
    wavl_node<k,v> *parent;
    do
    {
        if (*key == *node->key) return node;
        parent = node;
        node = *key < *node->key ? node->left : node->right;
    } while (node);

    wavl_node<k,v> *newNode = WavlTreeNewNode(arena, parent, key);
    ++tree->size;
    (*key < *parent->key ? parent->left : parent->right) = newNode;

    if (!parent->rank)
    {
        parent->rank = 1;
        WavlTreeBalance(tree, parent);
    }

    return newNode;
}

template
<typename k, typename v>
wavl_node<k,v> **
WavlTreeGetNodes(wavl_node<k,v> *node, wavl_node<k,v> **nodes)
{
    if (node) 
    { 
        nodes = WavlTreeGetNodes(node->left, nodes); 
        *(nodes++) = node;
        nodes = WavlTreeGetNodes(node->right, nodes); 
    }

    return nodes; 
}

template
<typename k, typename v>
wavl_node<k,v> **
WavlTreeGetNodes(wavl_tree<k,v> *tree, memory_arena *arena)
{
    wavl_node<k,v> **nodes = (wavl_node<k,v> **)PushSize_(arena, tree->size * sizeof(wavl_node<k,v> **));
    WavlTreeGetNodes(tree->root, nodes);
    return nodes;
}

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
    ll *refNames;
    wavl_tree<u64_string, basic_stats> *basicStats;
    wavl_tree<u64_string, wavl_tree<s32, wavl_tree<u64_string, ll>>> *moleculeData;
};

u08 LinkStats(link_stats_run_args *args, link_stats_return_data &data);

