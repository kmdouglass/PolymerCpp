/* Classes and functions pertaining to WormlikeChain generation
 * and analysis. */

#ifndef H_WORMLIKECHAIN
#define H_WORMLIKECHAIN

#include <iostream>       // Input/output
#include <Eigen/Core>     // Linear algebra
#include <Eigen/Geometry> // Cross product
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <random>         // Generating random numbers
#include <stdexcept>      // Throwing exceptions

#include "Misc.h"
#include "RgDict.h"
using namespace std;

const double pi = 3.1415926;


class Path
{   
public:
    int numPaths;
    vector<double> pathLength;
    double linDensity;
    double segConvFactor;
    Eigen::Vector3d initPoint;
    Eigen::Vector3d * points;
    std::vector<Eigen::Vector3d> path;

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
};

class WormlikeChain: public Path
/* A 3D wormlike chain.
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
 * locPrecision: double
 *      localization precision used for bumping the chain locations
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

    WormlikeChain(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_segConvFactor, double in_locPrecision, 
                  Eigen::Vector3d * in_initPoint);
    /* Constructor */

    void makePath(double in_pathLength);
    /*  Create the wormlike chain.
     *  
     *  The wormlike chain is created by first choosing the sizes of
     *  the small, random displacements in a plane tangent to a point
     *  on the surface of the unit sphere defined by the vector
     *  currPoint. The distribution for the sizes is given by the
     *  Boltzmann statistics for a semiflexible rod bending by a given
     *  angle due to interaction with its thermal environment.
     *  
     *  A random direction in this tangent plane is chosen by randomly
     *  and uniformly generating a vector on the unit sphere, taking
     *  its cross product with the currPoint vector, and normalizing
     *  the cross product. This cross product is multiplied by the
     *  size of the displacement found previously to generate the
     *  displacement vector.
     *  
     *  After displacing the currPoint vector into the tangent plane,
     *  the point in the plane is back projected onto the unit sphere
     *  to find the vector representing the next step in the polymer
     *  walk.
     *  
     *  This process is repeated until a number of vectors determined
     *  by numSegments representing a random walk on the surface of a
     *  sphere are generated. These vectors are cumulatively summed at
     *  the end to produce the final path field, which is the
     *  trajectory of the polymer.
     *
     *  Initial point to start polymer is determined by initPoint,
     *  which is set when initializing the class WormlikeChain. */

    void makeNewPath(double pathLength);
    /* Clears current path and makes a new one.
     * First point is determined by initPoint. */

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
};

RgDict * parSimChain(WormlikeChain * chain);
/*   Pimary processing for-loop to be parallelized.
  * 
  * parSimChain(data) is the most intensive part of the simulation. It
  * is a function applied to a WormlikeChain instance and repeatedly
  * calculates new conformations and gyration radii for those
  * conformations. Each WormlikeChain instance was defined with a
  * different persistence length.
  * 
  * Parameters
  * ----------
  * chainDict : WormlikeChainDicts
  *     A dictionary which contains information about
  *     one WormlikeChain instance
  * 
  * Returns
  * -------
  * RgDict : dictionary class
  *     Dictionary with Rg and RgBump keys containing the gyration
  *     radii for the chain and its sampled version.
  */

 double theoreticalWLCRg(double c, double Lp, double N);
 /* Return the theoretical value for the gyration radius.
  *
  *  Parameters
  *  ----------
  *  c : double
  *      The linear density of base pairs in the chain.
  *  Lp : double
  *      The persistence length of the wormlike chain.
  *  N : double
  *      The number of base pairs in the chain.
  *
  *  Returns
  *  -------
  *  meanRg : double 
  *     The mean gyration radius of a theoretical wormlike chain. */

#endif