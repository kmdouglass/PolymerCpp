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
 */
{
public:
	double defaultWeight; // always 1 in case of SAWLC
	SAWLC(double in_pathLength, 
               double in_persisLength,
               double in_linkDiameter, 
               Eigen::Vector3d * in_initPoint);
	/* Constructor */

	void makePath(double in_pathLength);
  // Simulates one instance of the path with length in_pathLength.

  virtual double getDefaultWeight();
  // SAWLC always has weight 1.

	virtual double getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where);
  // SAWLC always has weight 1.

	void getPossibleCollisions(vector<Eigen::Vector3d> & cumulativePath,
                    vector<int> & collisionPositions,
                    int where);
  /* Scans the already generated part of chain for possible collisions
   * with next step - calculates maximal extent of the excluded volume
   * sphere in next step and adds the index numbers of all points that
   * are closer than that extent to the current end of chain.       */

	bool checkCollision(Eigen::Vector3d * nextPoint,
                           vector<Eigen::Vector3d> & cumulativePath,
                           vector<int> & collisionPositions,
                           int where);
  /* Checks if the generated new step collides with some of the candidates
   * which were identified by getPossibleCollisions. */

	virtual double computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions,
                       int where);
  // SAWLC always has weight 1.
};


class SACollector : public Collector
{
public:
    vector<SAWLC*> myChains;
    // Constructor, initializes the variables, opens file for read/write
    // but doesn't start it.
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
    // Destructor, collects garbage.

    virtual SAWLC * getChainPointer(int i);
    // Returns pointer to the i-th chain in its index so it can be simulated
    // and its properties calculated and saved.
};


#endif