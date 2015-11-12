/* Collector class which gathers WormlikeChains */

#ifndef H_WLCCOLLECTOR
#define H_WLCCOLLECTOR

#include <iostream>       // Input/output
#include <Eigen/Core>     // Linear algebra
#include <Eigen/Geometry> // Cross product
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <stdexcept>      // Throwing exceptions
#include <fstream>        // Saving/Loading data
#include <limits>         // Limits of data types
#include <iomanip>        // std::setprecision

#include "WormlikeChain.h"
#include "RgDict.h"
#include "Misc.h"

class WLCCollector
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

    int numPaths;
    std::vector<double> pathLength;
    std::string nameDB;
    std::ofstream fileDB;
    double segConvFactor;
    std::vector<double> linDensity;
    std::vector<double> persisLength;
    double locPrecision;
    bool fullSpecParam;

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


    void startCollector();
    /* Simulates the WLCs and collects the information about them.
     * Then it saves the information into the specified database */
};

#endif