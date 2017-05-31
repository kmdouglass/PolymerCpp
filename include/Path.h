#ifndef H_PATH
#define H_PATH

#include <iostream>       // Input/output
#include <Eigen/Core>     // Linear algebra
#include <Eigen/Geometry> // Cross product
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <random>         // Generating random numbers
#include <stdexcept>      // Throwing exceptions
#include <iomanip>        // std::setprecision

#include "Misc.h"
#include "RgDict.h"
using namespace std;

const double pi = 3.1415926;


class Path
{   
public:
    double pathLength;
    double linkDiameter; // set to 0 if not SAWLC
    double persisLength;

    double chainWeight; // set to 1 if not Rosenbluth
    Eigen::Vector3d initPoint;
    Eigen::Vector3d * points;
    std::vector<Eigen::Vector3d> path;

    // constructor
    Path(double in_pathLength, 
        Eigen::Vector3d * in_initPoint);

    Eigen::Vector3d * randPointSphere(int numPoints);

	Eigen::Vector3d * randPointSphere();
    /* Randomly select points from the surface of a sphere.
     * Parameters
     * ----------
     * numPoints : int
     * The number of points to return.
     *
     * Returns
     * -------
     * points : pointer to array of Eigen::Vector3d of length numPoints
     *    The x, y, and z coordinates of each point on the sphere.
     *
     * References
     * ----------
     * [1] Weisstein, Eric W. "Sphere Point Picking." From
     * MathWorld--A Wolfram Web
     * Resource. http://mathworld.wolfram.com/SpherePointPicking.html
     */

    double computeRg();
    /* Compute the radius of gyration of a path.
     *
     * computeRg() calculates the radius of gyration of the WLC
     * object. The Rg is returned as a single number.
     *
     * Returns
     * -------
     * Rg : double
     *     The radius of gyration of the path object. */

    void bumpPoints(double locPrecision);
    /*   Bumps the points of this chain in a random direction in 3D.
     *   
     *   Parameters
     *   ----------
     *   locPrecision : double
     *       The localization precision of the measurement. This is the
     *       standard deviation of the Gaussian distribution
     *       determining the bump distances.
     */

    void makeNewPath(double pathLength);
    /* Clears current path and makes a new one.
     * First point is determined by initPoint. */

    virtual void makePath(double);

 /*TODO: needs to be restructured like the rest of the project.
  * 
  * RgDict * parSimChain();
  * Primary processing for-loop to be parallelized.
  * 
  * parSimChain(data) is the most intensive part of the simulation. It
  * is a function applied to a Path instance and repeatedly
  * calculates new conformations and gyration radii for those
  * conformations.
  * 
  * Returns
  * -------
  * RgDict : dictionary class
  *     Dictionary with Rg and RgBump keys containing the gyration
  *     radii for the chain and its sampled version.
  */
};



double theoreticalWLCRg(double Lp, double Lc);
 /* Return the theoretical value for the gyration radius.
  *
  * Parameters
  * ----------
  * Lp
  *     Persistence length of the chain.
  * Lc   
  *     Contour length of the chain.
  *
  *  Returns
  *  -------
  *  Rg2 
  *     The mean gyration radius (squared) of a theoretical wormlike chain. */


class Collector
{
public:
	int numPaths;
    std::vector<double> pathLength;
    std::string nameDB; // name of output file
    std::ofstream fileDB; // stream to output file
    double segConvFactor;
    std::vector<double> linDensity;
    std::vector<double> persisLength;
    std::vector<double> linkDiameter;
    double locPrecision;
    bool fullSpecParam;


    Collector(int in_numPaths,
              std::vector<double> & in_pathLength,
              std::vector<double> & in_linDensity,
              std::vector<double> & in_persisLength,
              std::string in_nameDB,
              double in_segConvFactor = 1.0,
              double in_locPrecision = 0.0,
              bool in_fullSpecParam = false);
    /* Initializes the collector - loads all values to member variables. */

    void startCollector();
    /* Runs the simulation and saves results to the "nameDB" file. */

    virtual Path * getChainPointer(int i);
    /* Used for passing different types of chain to the parallelized simulation
     * (prevents code duplicity)
     */
};


#endif