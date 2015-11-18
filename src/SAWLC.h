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
#include <fstream>        // Saving/Loading data
#include <limits>         // Limits of data types
#include <iomanip>        // std::setprecision

#include "Path.h"
#include "RgDict.h"
#include "Misc.h"

using namespace std;


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
	double defaultWeight;
	SAWLC(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_linkDiameter, double in_segConvFactor, 
                  double in_locPrecision, Eigen::Vector3d * in_initPoint);
	/* Constructor */

	void makePath(double in_pathLength);

	double getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where);

	void getPossibleCollisions(vector<Eigen::Vector3d> & cumulativePath,
                    vector<int> & collisionPositions,
                    int where);

	bool checkCollision(Eigen::Vector3d nextPoint,
                           vector<Eigen::Vector3d> & cumulativePath,
                           vector<int> & collisionPositions,
                           int where);

	double computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions,
                       int where);
};


class SACollector : public Collector
{
public:
    vector<SAWLC*> myChains;
    // Constructor, initializes the variables, opens file for read/write
    // and starts the collector
    SACollector(int in_numPaths,
                 std::vector<double> & in_pathLength,
                 std::vector<double> & in_linDensity,
                 std::vector<double> & in_persisLength,
                 std::vector<double> & in_linkDiameter,
                 std::string in_nameDB,
                 double in_segConvFactor = 1.0,
                 double in_locPrecision = 0.0,
                 bool in_fullSpecParam = false);

    ~SACollector();


    //void startCollector();
    /* Simulates the WLCs and collects the information about them.
     * Then it saves the information into the specified database */

    SAWLC * getChainPointer(int i);
};


#endif