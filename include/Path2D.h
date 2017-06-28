#ifndef H_PATH2D
#define H_PATH2D

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

class Path2D
{   
public:
    int pathLength;
    double linkDiameter; // set to 0 if not self-avoiding
    double persisLength;

    double chainWeight;  // set to 1 if not implementing Rosenbluth method
    Eigen::Vector2d initPoint;
    Eigen::Vector2d * points;
    std::vector<Eigen::Vector2d> path;

    // constructor
    Path2D(int in_pathLength, 
          Eigen::Vector2d * in_initPoint);

    double computeRg();
    /* Compute the radius of gyration of a 2D path.
     *
     * computeRg() calculates the radius of gyration of a 2D path.
     * The Rg is returned as a single number.
     *
     * Returns
     * -------
     * Rg : double
     *     The radius of gyration of the path object. */

    void makeNewPath(double pathLength);
    /* Clears current path and makes a new one.
     * First point is determined by initPoint. */

    virtual void makePath(double);

};

#endif