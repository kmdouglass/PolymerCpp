#ifndef H_PYUTILS
#define H_PYUTILS

#include <Python.h> // Must be first
#include <vector>
#include <stdexcept>

using namespace std;

PyObject* vectorToList_Float(const vector<double> &data);

PyObject* vectorToTuple_Float(const vector<double> &data);

vector<double> listTupleToVector_Float(PyObject* incoming);
#endif