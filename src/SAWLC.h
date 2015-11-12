#ifndef H_SAWLC
#define H_SAWLC

#include <iostream>       // Input/output
#include <Eigen/Core>     // Linear algebra
#include <Eigen/Geometry> // Cross product
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <vector>         // Vectors for storing data
#include <random>         // Generating random numbers
#include <stdexcept>      // Throwing exceptions

#include "WormlikeChain.h"

using namespace std;

extern double pi;

class SAWLC: public Path
/* A self-avoiding wormlike chain
 * 
 * Member variables:
 * -----------------
 * numPaths: int
 *     number of paths to be simulated
 * pathLength: vector<double> of length numPaths
 *     for each path contains its genomic length
 * linDensity: double
 *      linear density of the chain 
 * persisLength: double
 *     The persistence length in units of chain segments.
 * linkDiameter: double
 * 	   The diameter of a link
 * initPoint : pointer to Eigen::Vector3d
 *     The coordinates of the first point of the chain.
 * path : <vector> of Eigen::Vector3d points
 *     vector, whose elements are 3d vectors describing
 *     endpoints of the segments comprising the path
 */
{
public:
	double persisLength;
	double locPrecision;
	SAWLC(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_linkDiameter, double in_segConvFactor, 
                  double in_locPrecision, Eigen::Vector3d * in_initPoint);
	/* Constructor */
};




#endif