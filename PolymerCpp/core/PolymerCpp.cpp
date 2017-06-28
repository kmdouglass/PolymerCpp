#include "Python.h"

#include <iostream>       // Input/output
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <ctime>          // Timing
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <random>         // Generating random numbers
#include <chrono>         // Timing

using namespace std;

#include "WLC.h"
#include "WLC2D.h"
#include "SAWLC.h"
#include "SAWLC_Rosenbluth.h"
#include "PyUtils.h"

extern std::minstd_rand randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

static PyObject * getWLCrgs(PyObject *self, PyObject *args)
{
	seedRandom();

	//Parse input into parameters
        int    numPaths;
	double pathLength;
	double persisLength;
    if (!PyArg_ParseTuple(args, "idd", &pathLength, &persisLength))
        return NULL;

    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    // Initialize the chain
    WLC chain = WLC(pathLength, 
                    persisLength,
                    &startDir);

    // Initialize vectors for storing the gyration radii
    vector<double> results;

    // Generate all chains and store the results of calculation
    for (unsigned int i=0; i<numPaths; i++)
    {
    	chain.makePath(pathLength);
    	results.push_back(chain.computeRg());
    }

    
    return vectorToTuple_Float(results);
}

static PyObject * getWLC(PyObject *self, PyObject *args)
{
	seedRandom();

	//Parse input
	double pathLength;
	double persisLength;
    if (!PyArg_ParseTuple(args, "dd", &pathLength, &persisLength))
        return NULL;

    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    WLC chain = WLC(pathLength, 
                    persisLength,
                    &startDir);
    chain.makePath(pathLength);


    // Store all 3D vector coordinates in a single 1D vector (yikes)
    /* Format:
     * (1x, 1y, 1z, 2x, 2y, 2z, 3x, 3y, 3z, ...)
     */
    vector<double> chainPoints;
    for (unsigned int j=0; j<chain.path.size(); j++)
    {
    	for (unsigned int i=0; i <3 /*EPFL!*/; i++)
    	{
    		chainPoints.push_back(chain.path.at(j)[i]);
    	}
    }
    // Let Python deal with reshaping the array back.
    return vectorToTuple_Float(chainPoints);
}

static PyObject * getWLC2D(PyObject *self, PyObject *args)
{
    seedRandom();
    
    int pathLength;
    double persisLength;
    if (!PyArg_ParseTuple(args, "id", &pathLength, &persisLength)) {
        return NULL;
    }
    
    Eigen::Vector2d startDir(1.0, 0.0);
    
    WLC2D chain = WLC2D(pathLength,
                        persisLength,
                        &startDir);
    chain.makePath(pathLength);
    
    // Store all 2D vector coordinates in a single 1D vector
    /* Format:
     * (1x, 1y, 2x, 2y, 3x, 3y, ...)
     */
    vector<double> chainPoints;
    for (unsigned int j=0; j<chain.path.size(); j++)
    {
    	for (unsigned int i = 0; i < 2; i++)
    	{
    		chainPoints.push_back(chain.path.at(j)[i]);
    	}
    }
    // Let Python deal with reshaping the array back.
    return vectorToTuple_Float(chainPoints);
    
}

static PyObject * getSAWLCrgs(PyObject *self, PyObject *args)
// Same as getWLCrgs, just added another input parameter - linkDiameter, and
// changed the type of chain from WLC to SAWLC.
{
	seedRandom();
	//Parse input
        int numPaths;
	double pathLength;
	double persisLength;
	double linkDiameter;
    if (!PyArg_ParseTuple(args, "iddd", &pathLength, &persisLength, &linkDiameter))
        return NULL;

    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    SAWLC chain = SAWLC(pathLength, 
                        persisLength,
                        linkDiameter,
                        &startDir);
    vector<double> results;
    for (unsigned int i=0; i < numPaths; i++)
    {
    	chain.makePath(pathLength);
    	results.push_back(chain.computeRg());
    };

    return vectorToTuple_Float(results);
    
}


static PyObject * getSAWLC(PyObject *self, PyObject *args)
// Same as getWLC, just added another input parameter - linkDiameter, and
// changed the type of chain from WLC to SAWLC.
{
	seedRandom();
	//Parse input
	double pathLength;
	double persisLength;
	double linkDiameter;
    if (!PyArg_ParseTuple(args, "ddd", &pathLength, &persisLength, &linkDiameter))
        return NULL;

    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    SAWLC chain = SAWLC(pathLength, 
                        persisLength,
                        linkDiameter,
                        &startDir);
    chain.makePath(pathLength);

    vector<double> chainPoints;
    for (unsigned int j=0; j<chain.path.size(); j++)
    {
    	for (unsigned int i=0; i <3 /*EPFL!*/; i++)
    	{
    		chainPoints.push_back(chain.path.at(j)[i]);
    	}
    }

    return vectorToTuple_Float(chainPoints);
} 

static PyMethodDef PolymerCppCoreMethods[] = 
{
    {"getWLCrgs",  getWLCrgs, METH_VARARGS,
     "Generate an infinitesimally-thin wormlike chain in three dimensions."},
    {"getWLC",  getWLC, METH_VARARGS,
     "Get the whole chain of certain parameters."},
    {"getWLC2D", getWLC2D, METH_VARARGS,
     "Generate an infinitesimally-thin wormlike chain in two dimensions."}, 
    {"getSAWLCrgs",  getSAWLCrgs, METH_VARARGS,
     "Get multiple SAWLC radii of gyration of certain parameters."},
    {"getSAWLC",  getSAWLC, METH_VARARGS,
     "Get the whole self-avoiding chain of certain parameters."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef PolymerCppCore =
{
    PyModuleDef_HEAD_INIT,
    "PolymerCppCore",
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    PolymerCppCoreMethods
};

PyMODINIT_FUNC
PyInit_PolymerCppCore(void)
{
  return PyModule_Create(&PolymerCppCore);
}

int main(int argc, char *argv[])
{
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }

    /* Add a built-in module, before Py_Initialize */
    PyImport_AppendInittab("PolymerCppCore", PyInit_PolymerCppCore);

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(program);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Optionally import the module; alternatively,
       import can be deferred until the embedded script
       imports it. */
    PyImport_ImportModule("PolymerCppCore");

    PyMem_RawFree(program);
    return 0;
}
