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
#include <chrono>         // Timing
#include <stdexcept>      // Throwing exceptions
#include <fstream>        // Saving/Loading data
#include <limits>         // Limits of data types
#include <iomanip>        // std::setprecision


#define pi 3.1415926

using namespace std;

void seedRandom();


// defines engine and distribution for random number generation
unsigned seed; 
std::default_random_engine randGenerator(1u);
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
    int numPaths;
    vector<double> pathLength;
    double linDensity;
    double persisLength;
    double segConvFactor;
    double locPrecision;
    Eigen::Vector3d initPoint;
    std::vector<Eigen::Vector3d> path;

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

    double computeRg(int dimensions = 3);
    /* Compute the radius of gyration of a path.
     *
     * computeRg() calculates the radius of gyration of the WLC
     * object. The Rg is returned as a single number.
     *
     * Parameters
     * ----------
     * dimensions : int
     *     Compute the two dimensional or three dimensional radius of
     *     gyration (Default: dimensions = 3)
     *
     * Returns
     * -------
     * Rg : double
     *     The radius of gyration of the path object. */
};

/*class WormlikeChainDict
 * A container class for the WormlikeChain instance 
 * 
 * Member variables:
 * -----------------
 * chain: pointer to a new WormlikeChain instance
 *      WLC instance has to be destroyed manually
 *     
 * numSegments: vector<double>
 *      Number of segments to simulate for each chain iteration
 * locPrecision: double
 *      localization precision used for bumping the chain locations 
 * linDensity: double
 *      linear density of the chain
 * persisLength: double
 *      persistence length of the chain

{
public:
    
    WormlikeChain * chain;
    vector<double> numSegments;
    double locPrecision;
    double linDensity;
    double persisLength;
    double segConvFactor;

    // Constructor
    WormlikeChainDict(WormlikeChain * in_chain,
                      vector<double> * in_numSegments,
                      double in_locPrecision,
                      double in_linDensity,
                      double in_persisLength,
                      double in_segConvFactor);
    // Destructor, deletes the WLC instance
    ~WormlikeChainDict();
};
*/
class RgDict
/* A container for gathered data of gyration radii
 *
 * Member variables:
 * -----------------
 * Rg: <vector> of doubles
 *      Contains gyration radii for the chain
 * RgBump: <vector> of doubles
 *      Contains gyration radii for the sampled version 
 * linDensity: double
 *      linear density of this chain
 * persisLength: double
 *      persistence length of this chain
 * ALL ARE CONVERTED TO USER DEFINED UNITS UPON CREATION */
{
public:

    std::vector<double> Rg;
    std::vector<double> RgBump;
    double pathLength;
    double linDensity;
    double persisLength;
    double segConvFactor;

    // Constructor
    RgDict(std::vector<double> * in_Rg, std::vector<double> * in_RgBump,
        double in_pathLength, double in_linDensity, double in_persisLength,
        double in_segConvFactor, bool convert = false);

    void addToDBfile(std::ofstream & fileDB);
    /* Adds the RgDict to database file in following format:
     * linear density
     * persistence length
     * number of paths simulated (N)
     * N points of Rg doubles
     * N points of RgBump doubles */

    void readFromDBfile(std::ifstream & fileDB);
    /* Reads from DB file - NOT IMPLEMENTED YET */

};

vector<RgDict*> * readFromDBfile(std::ifstream & fileDB);
/* Reads from DB file.
 *
 * Returns a pointer to new vector of pointers to RgDict
 * instances, which are loaded with the data from file.
 *
 * Parameters:
 * -----------
 * fileDB: std::ifstream instance
 *      file opened with read privileges, previously created
 *      by this program
 *
 * Returns:
 * --------
 * data: pointer to new vector<RgDict*> */

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
     * outVector : vector of doubles
     *     Resulting segments
     * inVector : vector of doubles
     *     The parameters to convert into segments
     * segConvFactor: double
     *     Conversion factor betwee the user units and path segments
     * multiplyBool : bool
     *     Multiply or divide by the conversion factor
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