/* Classes and functions for storing data about gyration
 * radii, loading and saving them to files.
 */

#ifndef H_RGDICT
#define H_RGDICT

#include <iostream>       // Input/output
#include <cmath>          // General math operations
#include <cstdlib>        // Standard library
#include <string>         // Strings
#include <vector>         // Vectors for storing data
#include <fstream>        // Saving/Loading data
#include <stdexcept>      // Throwing exceptions
#include <sstream>

#include "Misc.h"

using namespace std;

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
    std::vector<double> Wt;
    std::vector<double> Rg;
    std::vector<double> RgBump;
    double pathLength;
    double linDensity;
    double persisLength;
    double linkDiameter;
    double segConvFactor;

    // Constructor
    RgDict(std::vector<double> & in_Rg, std::vector<double> & in_RgBump,
        std::vector<double> & in_Wt,
        double in_pathLength, double in_linDensity, double in_persisLength,
        double in_linkDiameter, double in_segConvFactor, bool convert = false);

    void addToDBfileFull(std::ofstream & fileDB);
    /* Adds the RgDict to database file in following format:
     * linear density
     * persistence length
     * link diameter
     * number of paths simulated (N)
     * N points of Rg doubles
     * N points of RgBump doubles */
     void addToDBfileShort(std::ofstream & fileDB);
    /* Adds the RgDict to database file in following format:
     * linear density
     * persistence length
     * link diameter
     * number of paths simulated (N)
     * mean of Rg
     * mean of RgBump */


};

#endif