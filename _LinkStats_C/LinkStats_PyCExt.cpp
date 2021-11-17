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
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
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

#define ArrayCount(array) (sizeof(array) / sizeof(array[0]))
#define ForLoop(n) for (u32 index = 0; index < (n); ++index)
#define ForLoop64(n) for (u64 index = 0; index < (n); ++index)
#define ForLoop2(n) for (u32 index2 = 0; index2 < (n); ++index2)
#define ForLoop3(n) for (u32 index3 = 0; index3 < (n); ++index3)
#define ForLoopN(i, n) for (u32 i = 0; i < (n); ++i)
#define TraverseLinkedList(startNode, type) for (type *(node) = (startNode); node; node = node->next)
#define TraverseLinkedList2(startNode, type) for (type *(node2) = (startNode); node2; node2 = node2->next)
#define TraverseLinkedList3(startNode, type) for (type *(node3) = (startNode); node3; node3 = node3->next)

#define ModuleName "_LinkStats"
#define ModuleDescription "LinkStats C Ext"

static
PyObject *
BasicStatsToTuple (basic_stats *stats)
{
    PyObject *insertSizes = PyTuple_New((Py_ssize_t)stats->insertSizes.count);
    u32 index = 0;
    TraverseLinkedList(stats->insertSizes.head, ll_node) PyTuple_SET_ITEM(insertSizes, (Py_ssize_t)index++, PyLong_FromUnsignedLongLong(node->value));

    return PyTuple_Pack(9, insertSizes, PyLong_FromUnsignedLongLong(stats->totalReadLength), PyLong_FromUnsignedLongLong(stats->totalAlignments), PyLong_FromUnsignedLongLong(stats->totalDup), PyLong_FromUnsignedLongLong(stats->totalQCF), PyLong_FromUnsignedLongLong(stats->totalUnM), PyLong_FromUnsignedLongLong(stats->totalNoMI), PyLong_FromUnsignedLongLong(stats->totalNoBX), PyLong_FromUnsignedLongLong(stats->totalZeroMQ));
}

static
PyObject *
MoleculeAlignmentsToTuple (const ll *list)
{
    PyObject *tuple = PyTuple_New((Py_ssize_t)list->count);
    {
	u32 index = 0;
	TraverseLinkedList(list->head, ll_node)
	{
	    alignment *al = (alignment *)node->ptr;
	    PyTuple_SET_ITEM(tuple, (Py_ssize_t)index++, PyTuple_Pack(6, PyBool_FromLong((s32)al->haveMI), PyLong_FromLong(al->qual), PyLong_FromUnsignedLongLong(al->referenceStart), PyLong_FromUnsignedLongLong(al->referenceEnd), PyLong_FromUnsignedLongLong(al->queryLength), PyLong_FromLong(al->mi)));
	}
    }
    return tuple;
}

static
PyObject *
MoleculeMap2ToTuple (const std::map<std::pair<std::string, s32>, ll *> &map2)
{
    PyObject *tuple = PyTuple_New((Py_ssize_t)map2.size());
    {
	u32 index = 0;
	for (const auto& [bxmi, list] : map2) PyTuple_SET_ITEM(tuple, (Py_ssize_t)index++, PyTuple_Pack(2, PyTuple_Pack(2, Py_BuildValue("s", bxmi.first.c_str()), PyLong_FromLong(bxmi.second)), MoleculeAlignmentsToTuple(list)));
    }
    return tuple;
}

static
PyObject *
MoleculeMap1ToTuple (const std::map<s32, std::map<std::pair<std::string, s32>, ll *>> &map1)
{
    PyObject *tuple = PyTuple_New((Py_ssize_t)map1.size());
    {
	u32 index = 0;
	for (const auto& [tid, map2] : map1) PyTuple_SET_ITEM(tuple, (Py_ssize_t)index++, PyTuple_Pack(2, PyLong_FromLong(tid), MoleculeMap2ToTuple(map2)));
    }
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
	write(errorFD, ErrorBuffer, n); \
	write(errorFD, "\n", 1); \
	PyErr_SetString(PyExc_Exception, ErrorBuffer); \
	return 0;\
    } while (0)

    s32 logFD;
    u32 nThreads;
    char *samFileName;
    char *fastaRefFileName;
    char *overrideName;
    char *fallbackName;
    u08 useMI;
    {
	PyObject *objParams;
	static const char *kwlist[] = {"params", 0};
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char **)kwlist, &objParams)) Error("_LinkStats_C No input 'params'");

	PyObject *objErrorFD, *objLogFD, *objNumThreads, *objSamFileName, *objFastaRefFileName, *objOverrideName, *objFallbackName, *objUseMI;

	{
	    struct
		inputParam
		{
		    PyObject **obj;
		    char *name;
		};

	    inputParam params[] = {
		{
		    &objLogFD,
		    (char *)"log"
		},
		{
		    &objErrorFD,
		    (char *)"error"
		},
		{
		    &objNumThreads,
		    (char *)"num_threads"
		},
		{
		    &objSamFileName,
		    (char *)"sam_file_name"
		},
		{
		    &objFastaRefFileName,
		    (char *)"fasta_file_name"
		},
		{
		    &objOverrideName,
		    (char *)"override_name"
		},
		{
		    &objFallbackName,
		    (char *)"fallback_name"
		},
		{
		    &objUseMI,
		    (char *)"use_mi"
		}
	    };

	    ForLoop(ArrayCount(params))
	    {
		inputParam *param = params + index;
		PyObject *objName = Py_BuildValue("s", param->name);
		*(param->obj) = PyObject_GetAttr(objParams, objName);
		Py_DECREF(objName);
		if (!(*(param->obj))) Error("_LinkStats_C No parm %s", param->name);
	    }
	}

	if ((logFD = (s32)PyObject_AsFileDescriptor(objLogFD)) < 0) Error("_LinkStats_C Invalid logFD '%d'", logFD);
	if ((errorFD = (s32)PyObject_AsFileDescriptor(objErrorFD)) < 0) Error("_LinkStats_C Invalid errorFD '%d'", errorFD);

	nThreads = (u32)PyLong_AsUnsignedLong(objNumThreads);
	Py_DECREF(objNumThreads);

	samFileName = (char *)PyUnicode_1BYTE_DATA(objSamFileName);
	fastaRefFileName = objFastaRefFileName == Py_None ? 0 : (char *)PyUnicode_1BYTE_DATA(objFastaRefFileName);
	overrideName = objOverrideName == Py_None ? 0 : (char *)PyUnicode_1BYTE_DATA(objOverrideName);
	fallbackName = (char *)PyUnicode_1BYTE_DATA(objFallbackName);

	useMI = PyObject_IsTrue(objUseMI) ? 1 : 0;
	Py_DECREF(objUseMI);
    }

    memory_arena workingSet;
    CreateMemoryArena(workingSet, MegaByte(512));

    link_stats_run_args runArgs;
    runArgs.logFD = logFD;
    runArgs.numThreads = nThreads;
    runArgs.samFileName = samFileName;
    runArgs.fastaReferenceFileName = fastaRefFileName;
    runArgs.overrideName = overrideName;
    runArgs.fallbackName = fallbackName;
    runArgs.useMI = useMI;
    runArgs.arena = &workingSet;

    link_stats_return_data data;
    u08 runState;
    Py_BEGIN_ALLOW_THREADS;
    runState = LinkStats(&runArgs, data);
    Py_END_ALLOW_THREADS;
    if (!runState) Error("_LinkStats_C Run Error");

    PyObject *genomeLength = PyLong_FromUnsignedLongLong(data.genomeLength);
    
    PyObject *refNames = PyTuple_New((Py_ssize_t)data.refNames.size());
    ForLoop((u32)data.refNames.size()) PyTuple_SET_ITEM(refNames, (Py_ssize_t)index, Py_BuildValue("s", data.refNames[index].c_str()));

    PyObject *basicStats = PyTuple_New((Py_ssize_t)data.basicStats.size());
    {
	u32 index = 0;
	for (const auto& [id, stats] : data.basicStats) PyTuple_SET_ITEM(basicStats, (Py_ssize_t)index++, PyTuple_Pack(2, Py_BuildValue("s", id.c_str()), BasicStatsToTuple(stats)));
    }

    PyObject *moleculeData = PyTuple_New((Py_ssize_t)data.moleculeData.size());
    {
	u32 index = 0;
	for (const auto& [id, map] : data.moleculeData) PyTuple_SET_ITEM(moleculeData, (Py_ssize_t)index++, PyTuple_Pack(2, Py_BuildValue("s", id.c_str()), MoleculeMap1ToTuple(map)));
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
