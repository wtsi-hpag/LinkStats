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

struct
u64_string
{
    s32 id;
    u32 length;
    u64 string;
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
};

void MakeCopy(memory_arena *arena, char_string *string, u64_string **out);

void MakeCopy(memory_arena *arena, s32 *id, s32 **out);

template
<typename k, typename v>
struct
wavl_node
{
    union
    {
        wavl_node<k,v> *parent;
        wavl_node<k,v> *next;
    };
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
u08
WavlTreeNeedToRotateLeftStrong(wavl_node<k,v> *node)
{
    return(((node->rank == 1) && !node->left) || (node->rank > (node->left->rank + 1)));
}

template
<typename k, typename v>
u08
WavlTreeNeedToRotateRightStrong(wavl_node<k,v> *node)
{
    return(((node->rank == 1) && !node->right) || (node->rank > (node->right->rank + 1)));
}

template
<typename k, typename v>
u08
WavlTreeNeedToRotateLeftWeak(wavl_node<k,v> *node)
{
    return((!node->left) || (node->rank > (node->left->rank + 1)));
}

template
<typename k, typename v>
u08
WavlTreeNeedToRotateRightWeak(wavl_node<k,v> *node)
{
    return((!node->right) || (node->rank > (node->right->rank + 1)));
}

template
<typename k, typename v>
void
WavlTreeRotateLeft(wavl_tree<k,v> *tree, wavl_node<k,v> *node)
{
    wavl_node<k,v> *right = node->right;
    node->right = right->left;
    if (right->left) right->left->parent = node;
    right->parent = node->parent;
    (node->parent ? (node->parent->left == node ? node->parent->left : node->parent->right) : tree->root) = right;
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
    if (left->right) left->right->parent = node;
    left->parent = node->parent;
    (node->parent ? (node->parent->right == node ? node->parent->right : node->parent->left) : tree->root) = left;
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
        else
        {
            if (WavlTreeNeedToRotateLeftStrong(parent))
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
}

template
<typename k, typename v, typename l>
wavl_node<k,v> *
WavlTreeFindNode(memory_arena *arena, wavl_tree<k,v> *tree, l *key)
{
    wavl_node<k,v> **node;
    if (!(*(node = &tree->root)))
    {
        ++tree->size;
        return (tree->root = WavlTreeNewNode(arena, (wavl_node<k,v> *)0, key));
    }
    
    wavl_node<k,v> *parent;
    do
    {
        if (*key == *(*node)->key) return *node;
        parent = *node;
        node = *key < *(*node)->key ? &(*node)->left : &(*node)->right;
    } while (*node);

    ++tree->size;
    auto result = *node = WavlTreeNewNode(arena, parent, key);

    if (!parent->rank)
    {
        parent->rank = 1;
        WavlTreeBalance(tree, parent);
    }

    return result;
}

template
<typename k, typename v>
wavl_node<k,v> *
WavlTreeFreezeToLL(wavl_node<k,v> *node, wavl_node<k,v> *prevNode, wavl_node<k,v> **head)
{
    if (node) 
    { 
        prevNode = WavlTreeFreezeToLL(node->left, prevNode, head); 
        (prevNode ? prevNode->next : *head) = node;
        node->next = 0;
        prevNode = WavlTreeFreezeToLL(node->right, node, (wavl_node<k,v> **)0); 
    }

    return prevNode;
}

template
<typename k, typename v>
wavl_node<k,v> *
WavlTreeFreezeToLL(wavl_tree<k,v> *tree)
{
    wavl_node<k,v> *head = 0;
    WavlTreeFreezeToLL(tree->root, (wavl_node<k,v> *)0, &head);
    
    return head; 
}

template
<typename t>
struct
ll_node
{
    ll_node *next;
    t data;
};

template
<typename t>
struct
ll
{
    ll_node<t> *head;
    ll_node<t> *tail;
    u64 count;
};

template
<typename t>
void
NewLL(memory_arena *arena, ll<t> **out)
{
    ll<t> *list = (ll<t> *)PushSize_(arena, sizeof(ll<t>));
    list->count = 0;
    list->head = 0;
    list->tail = (ll_node<t> *)list;
    *out = list;
}

template
<typename t>
void
LLAddValue(ll<t> *list, t data, memory_arena *arena)
{
    ll_node<t> *newNode = (ll_node<t> *)PushSize_(arena, sizeof(ll_node<t>)); 
    newNode->data = data; 
    newNode->next = 0; 
    list->tail = (list->tail->next = newNode); 
    ++list->count; 
}

#define InsertSize_Median_Estimate_Histogram_Size 2048
struct
insertsize_histogram
{
    u64 hist[InsertSize_Median_Estimate_Histogram_Size];
};

insertsize_histogram *NewInsertSizeHistogram(memory_arena *arena);
u64 EstimateMedian(insertsize_histogram *hist);

struct
bit_array
{
   u08 *bits;
   u64 size;
};

bit_array *CreateBitArray(u64 size, memory_arena *arena);
void FillBitArray(bit_array *array, u64 from, u64 to);
u08 IsBitSet(bit_array *array, u64 pos);
ll<u64> *FindGaps(bit_array *array, memory_arena *arena);
