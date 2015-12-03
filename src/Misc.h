/* Random number generation:
 * - randGenerator, randUniformReal and randNormalReal
 *    have to be declared as extern before use
 *
 * Conversion of segments between simulation and user-defined units.
 */

#ifndef H_MISC
#define H_MISC

#include <iostream>
#include <random>         // Generating random numbers
#include <chrono>         // Timing
#include <vector>
#include <cmath>

using namespace std;



void seedRandom();
void printRandom();
/* Miscellaneous functions to test random number generation. */



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


double sum(vector<double> & vect);
    /* Sums components of vector. */
#endif