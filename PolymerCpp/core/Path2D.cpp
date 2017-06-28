#include "Path2D.h"

extern std::minstd_rand randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

Path2D::Path2D(int in_pathLength, 
              Eigen::Vector2d * in_initPoint)
{
    pathLength   = in_pathLength;
    linkDiameter = 0.0; // default value, overridden in SAWLC
    chainWeight  = 1.0; // default value, overridden in Rosenbluth
    initPoint    = *in_initPoint;
    path.push_back(initPoint);
}

void Path2D::makeNewPath(double in_pathLength)
{
    pathLength = in_pathLength;
    /* The next block is for the case, when the path
     * runs into a blind spot - after a lot of unsuccessful
     * attempts to get out of the spot, function returns
     * which just runs a new simulation.
     */
    do
    {
        path.clear();                   //clears existing path
        path.push_back(initPoint);      //sets first point
        makePath(in_pathLength);
    } while (path.size() < in_pathLength);
}

double Path2D::computeRg()
{
    double Rg = 0.0;
    double secondMoments[2] = {0.0, 0.0};
    double mean; double variance;
    /* Radius of gyration is equal to sqrt of sums of variances
     * in each axis.
     */
    for (int i=0; i<2; i++) {
        mean = 0.0; variance = 0.0;
        for (auto & point : path ) { // iterate over all points in path
            mean += point(i);
        }
        mean = mean / path.size();

        for (auto & point : path ) {
            variance += (mean-point(i))*(mean-point(i));
        }
        variance = variance / path.size();
        secondMoments[i] = variance;
    }

    Rg = sqrt(secondMoments[0] + secondMoments[1]);
    return Rg;
}

void Path2D::makePath(double in_pathLength) {
	cout << "Ran the wrong makePath!"; 
    // This function is overridden by specific chain algorithms.
}