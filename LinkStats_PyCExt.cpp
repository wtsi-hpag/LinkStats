#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_8_API_VERSION
#include <Python.h>
#include "numpy/arrayobject.h"

#include <LinkStats.hpp>

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
Main (PyObject *self, PyObject *args, PyObject *kwargs)
{
   Py_RETURN_NONE; 
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
	PyObject *obj = PyModule_Create(&Module);
	import_array();
	return(obj);
}
