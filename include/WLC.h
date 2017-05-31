/* Classes and functions pertaining to WLC generation
 * and analysis. */

#ifndef H_WLC
#define H_WLC

#include <iostream>       // Input/output
#include <Eigen/Core>     // Linear algebra
#include <Eigen/Geometry> // Cross product
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <random>         // Generating random numbers
#include <stdexcept>      // Throwing exceptions
#include <fstream>        // Saving/Loading data
#include <limits>         // Limits of data types
#include <iomanip>        // std::setprecision

#include "Path.h"
#include "RgDict.h"
#include "Misc.h"

class WLC: public Path
/* A 3D wormlike chain.
 * 
 * Member variables:
 * -----------------
 * pathLength: vector<double> of length numPaths
 *     for each path contains its genomic length
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

    WLC(double in_pathLength,
         double in_persisLength, 
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
     *  which is set when initializing the class WLC. */

    double getTheoreticalRg();
};


class WLCCollector : public Collector
/*   Collector for the wormlike chain.
 *   Creates random walk paths and collects their statistics.
 * 
 *   A Collector generates a user-defined number of random walk paths
 *   with possibly different lengths by sending the walk parameters to
 *   a Path object. After the path has been generated, the statistics
 *   that describe the path are collected and binned into a histogram.
 *
 *   Member variables:
 *   -----------------
 *   numPaths : int
 *       The number of paths to collect before stopping the simulation
 *   pathLength : vector of doubles
 *       The length of each simulated path in genomic length
 *   linDensity : vector of doubles
 *       The number of base pairs per user-defined unit of length
 *   persisLength : vector of doubles
 *       The path's persistence length in user-defined units of length
 *   segConvFactor : double (optional)
 *       Conversion factor between the user units and path segments
 *       (Default is 1.0)
 *   locPrecision : double (optional)
 *       Standard deviation of the Gaussian defining the effective
 *       system PSF. (Default is 0, meaning no bumps are made)
 *   fullSpecParam : bool (optional)
 *       Do linDensity and persisLength define all the parameter-space
 *       points to simulate, or do they instead define the points in a
 *       grid to be generated with meshgrid? In the first case, the
 *       number of points to simulate is equal to the length of
 *       persisLength OR linDensity, whereas in the second case it's
 *       equal to the number of points in persisLength TIMES the number
 *       of points in linDensity. (Default is false; the points will
 *       define a grid in this case).
 */
{
public:
    vector<WLC*> myChains;

    // Constructor, initializes the variables, opens file for read/write
    // and starts the collector
    WLCCollector(int in_numPaths,
                 std::vector<double> & in_pathLength,
                 std::vector<double> & in_linDensity,
                 std::vector<double> & in_persisLength,
                 std::string in_nameDB,
                 double in_segConvFactor = 1.0,
                 double in_locPrecision = 0.0,
                 bool in_fullSpecParam = false);

    ~WLCCollector();

    WLC * getChainPointer(int i);
    //void startCollector();
    /* Simulates the WLCs and collects the information about them.
     * Then it saves the information into the specified database */
};

#endif