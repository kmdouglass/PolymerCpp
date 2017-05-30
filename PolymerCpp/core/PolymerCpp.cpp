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
	int numPaths;
	double pathLength;
	double linDensity;
	double persisLength;
	double segConvFactor;
	int alsoBumped;
	double locPrecision;
    if (!PyArg_ParseTuple(args, "iddddid", &numPaths, &pathLength,
    	&linDensity, &persisLength, &segConvFactor, &alsoBumped,
    	&locPrecision))
        return NULL;

    /* Initialize pathLength vector. All paths have equal lengths because
     * I didn't figure out how to parse a Python tuple of variable length
     * into a C++ function.) */
	vector<double> pathLengths(numPaths, pathLength);
    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    // Initialize the chain
    WLC chain = WLC(numPaths, pathLengths, 
                  linDensity, persisLength,
                  segConvFactor, locPrecision, 
                  &startDir);

    // Initialize vectors for storing the gyration radii
    vector<double> results;
    vector<double> resultsBumped;

    // Generate all chains and store the results of calculation
    for (unsigned int i=0; i<numPaths; i++)
    {
    	chain.makePath(pathLengths.at(i));
    	results.push_back(chain.computeRg());
    	if (alsoBumped)
    	{
    		chain.bumpPoints(locPrecision);
    		resultsBumped.push_back(chain.computeRg());
    	}
    }

    /* If we also bumped the chain, concatenate the two results vectors,
	 * and add a delimiting negative number between them. This should be
	 * the only negative number in the results.
	 */
    if (alsoBumped)
    {
    	vector<double> resultsCombined;
    	resultsCombined.reserve( results.size() + resultsBumped.size() + 1u); // preallocate memory
		resultsCombined.insert( resultsCombined.end(), results.begin(), results.end() );
		resultsCombined.push_back(-1.0); // delimiting negative number
		resultsCombined.insert( resultsCombined.end(), resultsBumped.begin(), resultsBumped.end() );
		return vectorToTuple_Float(resultsCombined);
    }

    // Else return just the un-bumped results.
    else 
    {
    	return vectorToTuple_Float(results);
    }
}

static PyObject * getWLC(PyObject *self, PyObject *args)
{
	seedRandom();

	//Parse input
	double pathLength;
	double linDensity;
	double persisLength;
	double segConvFactor;
	int bumped;
	double locPrecision;
    if (!PyArg_ParseTuple(args, "ddddid", &pathLength,
    	&linDensity, &persisLength, &segConvFactor, &bumped,
    	&locPrecision))
        return NULL;

    // Same as previous function
	vector<double> pathLengths(1, pathLength);
    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    WLC chain = WLC(1, pathLengths, 
                  linDensity, persisLength,
                  segConvFactor, locPrecision, 
                  &startDir);
    chain.makePath(pathLength);


    if (bumped)
    {
    	chain.bumpPoints(locPrecision);
    }

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



static PyObject * getSAWLCrgs(PyObject *self, PyObject *args)
// Same as getWLCrgs, just added another input parameter - linkDiameter, and
// changed the type of chain from WLC to SAWLC.
{
	seedRandom();
	//Parse input
	int numPaths;
	double pathLength;
	double linDensity;
	double persisLength;
	double segConvFactor;
	int alsoBumped;
	double locPrecision;
	double linkDiameter;
    if (!PyArg_ParseTuple(args, "iddddidd", &numPaths, &pathLength,
    	&linDensity, &persisLength, &segConvFactor, &alsoBumped,
    	&locPrecision, &linkDiameter))
        return NULL;

	vector<double> pathLengths(numPaths, pathLength);
    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    SAWLC chain = SAWLC(numPaths, pathLengths, 
                  linDensity, persisLength,
                  linkDiameter,
                  segConvFactor, locPrecision, 
                  &startDir);
    vector<double> results;
    vector<double> resultsBumped;
    for (unsigned int i=0; i<numPaths; i++)
    {
    	chain.makePath(pathLengths.at(i));
    	results.push_back(chain.computeRg());
    	if (alsoBumped)
    	{
    		chain.bumpPoints(locPrecision);
    		resultsBumped.push_back(chain.computeRg());
    	}
    };

    if (alsoBumped)
    {
    	vector<double> resultsCombined;
    	resultsCombined.reserve( results.size() + resultsBumped.size() + 1u); // preallocate memory
		resultsCombined.insert( resultsCombined.end(), results.begin(), results.end() );
		resultsCombined.push_back(-1.0); // delimiting negative number
		resultsCombined.insert( resultsCombined.end(), resultsBumped.begin(), resultsBumped.end() );
		return vectorToTuple_Float(resultsCombined);
    }
    else 
    {
    	return vectorToTuple_Float(results);
    }
}


static PyObject * getSAWLC(PyObject *self, PyObject *args)
// Same as getWLC, just added another input parameter - linkDiameter, and
// changed the type of chain from WLC to SAWLC.
{
	seedRandom();
	//Parse input
	double pathLength;
	double linDensity;
	double persisLength;
	double segConvFactor;
	int bumped;
	double locPrecision;
	double linkDiameter;
    if (!PyArg_ParseTuple(args, "ddddidd", &pathLength,
    	&linDensity, &persisLength, &segConvFactor, &bumped,
    	&locPrecision, &linkDiameter))
        return NULL;

	vector<double> pathLengths(1, pathLength);
    Eigen::Vector3d startDir(1.0,0.0,0.0);
    
    SAWLC chain = SAWLC(1, pathLengths, 
                  linDensity, persisLength,
                  linkDiameter,
                  segConvFactor, locPrecision, 
                  &startDir);
    chain.makePath(pathLength);

    if (bumped)
    {
    	chain.bumpPoints(locPrecision);
    }

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

static PyMethodDef PolymerCppMethods[] = 
{
    {"getWLCrgs",  getWLCrgs, METH_VARARGS,
     "Get multiple WLC radii of gyration of certain parameters."},
    {"getWLC",  getWLC, METH_VARARGS,
     "Get the whole chain of certain parameters."},
     {"getSAWLCrgs",  getSAWLCrgs, METH_VARARGS,
     "Get multiple SAWLC radii of gyration of certain parameters."},
    {"getSAWLC",  getSAWLC, METH_VARARGS,
     "Get the whole self-avoiding chain of certain parameters."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef PolymerCpp =
{
    PyModuleDef_HEAD_INIT,
    "PolymerCpp",
    "",          /* module documentation, may be NULL */
    -1,          /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    PolymerCppMethods
};

PyMODINIT_FUNC
PyInit_PolymerCpp(void)
{
  return PyModule_Create(&PolymerCpp);
}

int main(int argc, char *argv[])
{
    wchar_t *program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL) {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        exit(1);
    }

    /* Add a built-in module, before Py_Initialize */
    PyImport_AppendInittab("PolymerCpp", PyInit_PolymerCpp);

    /* Pass argv[0] to the Python interpreter */
    Py_SetProgramName(program);

    /* Initialize the Python interpreter.  Required. */
    Py_Initialize();

    /* Optionally import the module; alternatively,
       import can be deferred until the embedded script
       imports it. */
    PyImport_ImportModule("PolymerCpp");

    PyMem_RawFree(program);
    return 0;
}
