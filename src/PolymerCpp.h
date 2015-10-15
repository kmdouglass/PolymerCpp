#ifndef H_POLYMERCPP
#define H_POLYMERCPP

#include <iostream>       // Input/output
#include <Eigen/Core>     // Linear algebra
#include <Eigen/Geometry> // Cross product
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <ctime>          // Timing
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <random>         // Generating random numbers
#include <stdexcept>      // Throwing exceptions

#define pi 3.1415926

using namespace std;

// defines engine and distribution for random number generation
std::default_random_engine randGenerator;
std::uniform_real_distribution<double> randUniformReal(0.0,1.0);
std::normal_distribution<double> randNormalReal(0.0,1.0);

class Path
{   
public:
    Eigen::Vector3d * points;
    Eigen::Vector3d * randPointSphere(int numPoints);
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
    void checkPath();
    /* Check if path has a correct number of dimensions. */
};

class WormlikeChain: public Path
/* A 3D wormlike chain.
 * 
 * Member variables:
 * -----------------
 * numSegments : float
 *     The number of segments in the chain. Must be >=1
 *     and need not be an integer.
 * pLength: int
 *     The persistence length in units of chain segments.
 * initPoint : pointer to Eigen::Vector3d
 *     The coordinates of the first point of the chain.
 * path : <vector> of Eigen::Vector3d points
 *     vector, whose elements are 3d vectors describing
 *     endpoints of the segments comprising the path
 */
{
public:

    double numSegments;
    int pLength;
    Eigen::Vector3d initPoint;
    std::vector<Eigen::Vector3d> path;

    WormlikeChain(double in_numSegments, int in_pLength,
                  Eigen::Vector3d * in_initPoint);
    /* Constructor */

    void makePath();
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

    void makeNewPath();
    /* Clears current path and makes a new one.
     * First point is determined by initPoint. */
};

class WormlikeChainDict
/* A container class for the WormlikeChain instance 
 * 
 * Member variables:
 * -----------------
 * chain: pointer to a new WormlikeChain instance
 *      WLC instance has to be destroyed manually
 *     
 * numSegments: double SHOULD BE VECTOR OF DOUBLE
 *      Number of segments to simulate for each chain iteration
 * locPrecision: double SHOULD ALSO BE VECTOR OF DOUBLES?
 *      localization precision used for bumping the chain locations */
{
public:
    
    WormlikeChain * chain;
    double numSegments;
    double locPrecision;

    // Constructor
    WormlikeChainDict(WormlikeChain * in_chain,
                      double in_numSegments,
                      double in_locPrecision);
    // Destructor, deletes the WLC instance
    ~WormlikeChainDict();
};

class RgDict
/* A container for gathered data of gyration radii
 *
 * Member variables:
 * -----------------
 * Rg: <vector> of doubles
 *      Contains gyration radii for the chain
 * RgBump: <vector> of doubles
 *      Contains gyration radii for the sampled version */
{
public:

    std::vector<double> Rg;
    std::vector<double> RgBump;

    // Constructor
    RgDict(std::vector<double> * in_Rg, std::vector<double> * in_RgBump);
};

RgDict * parSimChain(WormlikeChainDict * chainDict);
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

class Collector
/* Creates random walk paths and collects their statistics.
 * 
 * A Collector generates a user-defined number of random walk paths
 * with possibly different lengths by sending the walk parameters to
 * a Path object. After the path has been generated, the statistics
 * that describe the path are collected and binned into a histogram.
 * 
 * Member variables:
 * -----------------
 * numPaths: int
 *      Number of paths to collect before stopping the simulation
 * pathLength : vector of floats
 *     The length of each simulated path
 * nameDB: string
 *     Name of database
 * segConvFactor : float
 *     Conversion factor between the user units and path segments
 *     (Default is 1.0) */
{
public:
    int numPaths;
    std::vector<double> pathLength;
    std::string nameDB;
    double segConvFactor;

    // Constructor
    Collector();

    // Constructor, not used
    Collector(std::vector<double> & in_pathLength, 
          std::string in_nameDB,
          double in_segConvFactor = 1.0);
protected:
    void convSegments(std::vector<double> & outVector,
                      std::vector<double> & inVector,
                      double convFactor,
                      bool multiplyBool);
    /* Convert path parameters into segments.
     * Not the most efficient way,
     * but this is not a demanding operation.
     * 
     * Parameters
     * ----------
     * outVector : vector of floats
     *     Resulting segments
     * inVector : vector of floats
     *     The parameters to convert into segments
     * segConvFactor: float
     *     Conversion factor betwee the user units and path segments
     * multiplyBool : bool
     *     Multiply or divide by the conversion factor
     */
};

class WLCCollector: public Collector
/*   Collector for the wormlike chain.
 *
 *   Member variables:
 *   -----------------
 *   numPaths : int
 *       The number of paths to collect before stopping the simulation
 *   pathLength : vector of floats
 *       The length of each simulated path in genomic length
 *   linDensity : vector of floats
 *       The number of base pairs per user-defined unit of length
 *   persisLength : vector of floats
 *       The path's persistence length in user-defined units of length
 *   segConvFactor : float (optional)
 *       Conversion factor between the user units and path segments
 *       (Default is 1)
 *   locPrecision : float (optional)
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
    // Constructor, initializes the variables and starts the collector
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

    std::vector<double> linDensity;
    std::vector<double> persisLength;
    double locPrecision;
    bool fullSpecParam;
};

double computeRg(WormlikeChain * chain, int dimensions = 3);
/* Compute the radius of gyration of a path.
 *
 * computeRg() calculates the radius of gyration of a Path
 * object. The Rg is returned as a single number.
 *
 * Parameters
 * ----------
 * chain : Pointer to WormlikeChain instance
 *     The path is located at chain->path. The radius of gyration is computed
 *     from the endpoints of its individual segments.
 * dimensions : int
 *     Compute the two dimensional or three dimensional radius of
 *     gyration (Default: dimensions = 3)
 *
 * Returns
 * -------
 * Rg : double
 *     The radius of gyration of the path object. */

void bumpPoints(WormlikeChain * chain, double locPrecision);
/*Bumps the points in a random direction in 3D.
 *   
 *   Parameters
 *   ----------
 *   chain: pointer to WormlikeChain instance to be bumped
 *       This instance of WormlikeChain will have its points bumped
 *   locPrecision : double
 *       The localization precision of the measurement. This is the
 *       standard deviation of the Gaussian distribution
 *       determining the bump distances.
 */



class DateClass
/* Class which stores date and time in formatted strings at the time
 * of its initialization */
{
private:
    std::string m_date;
    std::string m_time;
    time_t tick;
    struct tm * timeinfo;

public:
    DateClass();
    std::string GetDate();
    std::string GetTime();
};
#endif