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

vector<RgDict*> * loadModel(vector<string> & dbNameList);
/* Loads the RgDicts from files.
 *
 * Parameters
 * ----------
 * dbNameList: vector of string
 *      Name(s) of the text files that contain the data
 * 
 * Returns
 * -------
 * data : pointer to vector of pointers to new RgDicts
 */

#endif