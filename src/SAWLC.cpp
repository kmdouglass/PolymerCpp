#include "SAWLC.h"

using namespace std;

extern std::default_random_engine randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

SAWLC::SAWLC(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_linkDiameter, double in_segConvFactor, 
                  double in_locPrecision, Eigen::Vector3d * in_initPoint)
{
	if (in_linDensity < 1.0/in_linkDiameter)
    {
        std::stringstream buffer;
        buffer << "The link diameter is smaller than the distance between links. "
            << std::endl << "Link diameter: " << in_linkDiameter << std::endl
            << "Distance between links: " << 1.0/in_linDensity << std::endl;
        throw std::out_of_range(buffer.str());
    }

    numPaths = in_numPaths;
    pathLength = in_pathLength;
    locPrecision = in_locPrecision;
    linDensity = in_linDensity;
    persisLength = in_persisLength;
    linkDiameter = in_linkDiameter;
    segConvFactor = in_segConvFactor;

    initPoint = *in_initPoint;
    path.push_back(initPoint);
}

void SAWLC::makePath(double in_pathLength)
{
    double numSegments = in_pathLength / linDensity;
    // check if numSegments is in valid range
    if (numSegments < 1)
    {
        std::stringstream buffer;
        buffer << "The number of segments must be greater than 1, but a "
            << "value of " << numSegments << "was supplied." << std::endl;
        throw std::out_of_range(buffer.str());
    }

	// split numSegments into integer and leftover
    int numSegInt; double numSegFrac;
    numSegInt = (int)numSegments;
    numSegFrac = numSegments - (double)numSegInt;
    //cout << "numSegFrac = " << numSegFrac << endl;
    // If numSegFrac is not zero, allocates memory also for leftover vector
    int numSeg = (abs(numSegFrac)>0.001) ? numSegInt+1 : numSegInt;

    double sigma = pow(2.0 / persisLength, 0.5);

    // Primary iterative loop for creating the chain
    Eigen::Vector3d currPoint = initPoint;
    Eigen::Vector3d * workingPath = new Eigen::Vector3d[numSeg];
    vector<double> stepWeights(numSeg) = 0.0;
    Eigen::Vector3d dispVector, nextPoint;
    double projDistance;
    workingPath[0] = currPoint;
    for (int i=1; i<numSeg; i++)
    {
    	// compute weight of step by integration over space
    		//check for possible collisions from other links
    			//get distance from each link
    				//if more than linkDiameter, OK

    				//if one less than linkDiameter, get interfering phase space,
    				//integrate it and subtract from all space

    				//if two less, check if their spaces conflict
    				//integrate them and subtract from all space

    				//if more than two, get available space and integrate over that 

    	// keep generating vectors until you get one which is ok

    	// assign it the weight

    	// move to next step

    }
}