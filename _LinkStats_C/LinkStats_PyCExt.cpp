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

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <LinkStats.hpp>

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wreserved-id-macro"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wcast-align"
#pragma GCC diagnostic ignored "-Wextra-semi-stmt"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wconditional-uninitialized"
#pragma GCC diagnostic ignored "-Wdouble-promotion"
#pragma GCC diagnostic ignored "-Wpadded"
#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#include "stb_sprintf.h"
#pragma clang diagnostic pop

#define ModuleName "_LinkStats"
#define ModuleDescription "LinkStats C Ext"

static
PyObject *
BasicStatsToTuple (basic_stats *stats)
{
    PyObject *insertSizes = PyTuple_New((Py_ssize_t)stats->insertSizes.count);
    TraverseLinkedList(stats->insertSizes.head) PyTuple_SET_ITEM(insertSizes, (Py_ssize_t)index, PyLong_FromUnsignedLongLong(node->data));

    return PyTuple_Pack(9, insertSizes, PyLong_FromUnsignedLongLong(stats->totalReadLength), PyLong_FromUnsignedLongLong(stats->totalAlignments), PyLong_FromUnsignedLongLong(stats->totalDup), PyLong_FromUnsignedLongLong(stats->totalQCF), PyLong_FromUnsignedLongLong(stats->totalUnM), PyLong_FromUnsignedLongLong(stats->totalNoMI), PyLong_FromUnsignedLongLong(stats->totalNoBX), PyLong_FromUnsignedLongLong(stats->totalZeroMQ));
}

static
PyObject *
MoleculeAlignmentsToTuple (ll<alignment *> *list)
{
    PyObject *tuple = PyTuple_New((Py_ssize_t)list->count);
    TraverseLinkedList(list->head) PyTuple_SET_ITEM(tuple, (Py_ssize_t)index, PyTuple_Pack(6, PyBool_FromLong((s32)node->data->haveMI), PyLong_FromLong(node->data->qual), PyLong_FromUnsignedLongLong(node->data->referenceStart), PyLong_FromUnsignedLongLong(node->data->referenceEnd), PyLong_FromUnsignedLongLong(node->data->queryLength), PyLong_FromLong(node->data->mi)));

    return tuple;
}

static
PyObject *
MoleculeMap2ToTuple (wavl_tree<u64_string, ll<alignment *>> *map2)
{
    PyObject *tuple = PyTuple_New((Py_ssize_t)map2->size);
    TraverseLinkedList(WavlTreeFreezeToLL(map2)) PyTuple_SET_ITEM(tuple, (Py_ssize_t)index, PyTuple_Pack(2, PyTuple_Pack(2, Py_BuildValue("s", (const char *)charU64String(node->key)), PyLong_FromLong(node->key->id)), MoleculeAlignmentsToTuple(node->value)));
    
    return tuple;
}

static
PyObject *
MoleculeMap1ToTuple (wavl_tree<s32, wavl_tree<u64_string, ll<alignment *>>> *map1)
{
    PyObject *tuple = PyTuple_New((Py_ssize_t)map1->size);
    TraverseLinkedList(WavlTreeFreezeToLL(map1)) PyTuple_SET_ITEM(tuple, (Py_ssize_t)index, PyTuple_Pack(2, PyLong_FromLong(*(node->key)), MoleculeMap2ToTuple(node->value)));

    return tuple;
}

static
PyObject *
Main (PyObject *self, PyObject *args, PyObject *kwargs)
{
    char ErrorBuffer[1024];
    s32 errorFD = -1;
    #define Error(message, ...) \
    do \
    { \
	u64 n = (u64)stbsp_snprintf(ErrorBuffer, sizeof(ErrorBuffer), message, ##__VA_ARGS__); \
	(void)write(errorFD, ErrorBuffer, n); \
	(void)write(errorFD, "\n", 1); \
	PyErr_SetString(PyExc_Exception, ErrorBuffer); \
	return 0;\
    } while (0)

    memory_arena workingSet;
    link_stats_run_args runArgs;
    {
	s32 logFD = -1;
	s32 nThreads = 1;
	char *samFileName = 0;
	char *fastaRefFileName = 0;
	char *overrideName = 0;
	char *fallbackName = 0;
	s32 useMI = 0;
	
	static const char *kwlist[] = {	"log",
					"error",
					"num_threads",
					"sam_file_name",
					"fasta_file_name",
					"override_name",
					"fallback_name",
					"use_mi",
					0};
	
	auto GetFD = (s32 (*)(PyObject *, void *))([](PyObject *obj, void *ptr)->s32
	{
	    s32 fd;
	    if ((fd = PyObject_AsFileDescriptor(obj)) < 0) return 0;
	    *((s32 *)ptr) = fd;
	    return 1;
	});
	
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|$O&O&iszzsp", (char **)kwlist, GetFD, &logFD,
											GetFD, &errorFD,
											&nThreads,
											&samFileName,
											&fastaRefFileName,
											&overrideName,
											&fallbackName,
											&useMI)) return 0;

	if (nThreads <= 0) Error("num_threads (%d) <= 0", nThreads);
	if (!samFileName) Error("sam_file_name not set");
	if (!fallbackName) Error("fallback_name not set");

	CreateMemoryArena(workingSet, MegaByte(512));

	runArgs.logFD = logFD;
	runArgs.numThreads = (u32)nThreads;
	runArgs.samFileName = PushU64String(samFileName, &workingSet);
	runArgs.fastaReferenceFileName = PushU64String(fastaRefFileName, &workingSet);
	runArgs.overrideName = PushU64String(overrideName, &workingSet);
	runArgs.fallbackName = PushU64String(fallbackName, &workingSet);
	runArgs.useMI = (u08)useMI;
	runArgs.arena = &workingSet;
    }
    
    link_stats_return_data *data;
    Py_BEGIN_ALLOW_THREADS;
    data = LinkStats(&runArgs);
    Py_END_ALLOW_THREADS;
    if (!data)
    {
	FreeMemoryArena(workingSet);
	Error("_LinkStats_C Run Error");
    }
    
    PyObject *genomeLength = PyLong_FromUnsignedLongLong(data->genomeLength);
    
    PyObject *refNames = PyTuple_New((Py_ssize_t)data->refNames->count);
    {
	TraverseLinkedList(data->refNames->head) PyTuple_SET_ITEM(refNames, (Py_ssize_t)index, Py_BuildValue("s", (const char *)charU64String(node->data)));
    }
    
    PyObject *basicStats = PyTuple_New((Py_ssize_t)data->basicStats->size);
    {
	TraverseLinkedList(WavlTreeFreezeToLL(data->basicStats)) PyTuple_SET_ITEM(basicStats, (Py_ssize_t)index, PyTuple_Pack(2, Py_BuildValue("s", (const char *)charU64String(node->key)), BasicStatsToTuple(node->value)));
    }
    
    PyObject *moleculeData = PyTuple_New((Py_ssize_t)data->moleculeData->size);
    {
	TraverseLinkedList(WavlTreeFreezeToLL(data->moleculeData)) PyTuple_SET_ITEM(moleculeData, (Py_ssize_t)index, PyTuple_Pack(2, Py_BuildValue("s", (const char *)charU64String(node->key)), MoleculeMap1ToTuple(node->value)));
    }

    FreeMemoryArena(workingSet);
    return PyTuple_Pack(4, genomeLength, refNames, basicStats, moleculeData); 
}

static
PyMethodDef
Methods[] =
{
    {ModuleName, (PyCFunction) Main, METH_VARARGS | METH_KEYWORDS, ModuleDescription},
    {NULL, NULL, 0, NULL}
};

static
struct
PyModuleDef
Module =
{
    PyModuleDef_HEAD_INIT,
    ModuleName,
    ModuleDescription,
    -1,
    Methods
};

PyMODINIT_FUNC
PyInit__LinkStats_C()
{
    return PyModule_Create(&Module);
}
