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
    int numSeg;
    numSeg = (int)numSegments;

    double sigma = pow(2.0 / persisLength, 0.5);

    // Primary iterative loop for creating the chain
    Eigen::Vector3d currPoint = initPoint;
    vector<Eigen::Vector3d> workingPath, cumulativePath;
    vector<double> stepWeights(numSeg,0.0);
    Eigen::Vector3d dispVector, nextPoint;
    double projDistance;
    workingPath.push_back(currPoint);
    cumulativePath.push_back(currPoint);
    vector<int> collisionPositions;
    double sigma = pow(2.0 / persisLength, 0.5);
    double angDisp, tanPlaneDisp;
    Eigen::Vector3d * randVec = NULL;
    for (int i=1; i<numSeg; i++)
    {
    	// compute weight of step by integration over space
        getPossibleCollisions(cumulativePath, collisionPositions);

        stepWeights[i] = getNextStepWeight(cumulativePath, collisionPositions); 

        getNextStep(workingPath, collisionPositions);
    	// keep generating vectors until you get one which is ok
        
        while (true) 
        {
            angDisp = sigma * randNormalReal(randGenerator);
            tanPlaneDisp = sin(angDisp);
            randVec = randPointSphere();
            dispVector = randVec.cross(currPoint);
            dispVector /= sqrt(dispVector.squaredNorm());
            dispVector = dispVector * tanPlaneDisp;
            projDistance = 1 - cos(angDisp);
            nextPoint = ((1 - projDistance) * currPoint) + dispVector;
            if (!checkCollision(nextPoint, cumulativePath, collisionPositions))
            {
                workingPath.push_back(nextPoint);
                cumulativePath.push_back(nextPoint + cumulativePath.back());
                currPoint = nextPoint;
                delete randVec;
                break;
            }
            delete randVec;
        }
    	// move to next step

    }
}

double SAWLC::getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions)
{
    // check for possible collisions from other links

    if (collisionPositions.size() == 0)
        return defaultWeight;
    else
        return computeIntegral(cumulativePath, collisionPositions);
}

void SAWLC::getPossibleCollisions(
                    vector<Eigen::Vector3d> & cumulativePath,
                    vector<int> & collisionPositions)
{
    collisionPositions.empty();
    Eigen::Vector3d distanceVector; double distance;
    Eigen::Vector3d endPoint = cumulativePath.back();
    for (int i = 0; i < cumulativePath.size()-2; i++)
    // not iterating over last two points, because they are always
    // colliding
    {
        distanceVector = endPoint - cumulativePath[i];
        distance = sqrt(distanceVector.squaredNorm());
        if (distance < (1.0 + linkDiameter))
        {
            collisionPositions.push_back(i);
        }
    }
}

bool SAWLC::checkCollision(Eigen::Vector3d nextPoint,
                           vector<Eigen::Vector3d> & cumulativePath,
                           vector<int> & collisionPositions)
{
    Eigen::Vector3d endPoint = nextPoint + cumulativePath.back();
    Eigen::Vector3d distanceVector; double distance;

    collisionPositions.push_back(cumulativePath.size() - 2);
    for (auto & position : collisionPositions)
    {
        distanceVector = endPoint - cumulativePath[position];
        distance = sqrt(distanceVector.squaredNorm());
        if (distance < (1.0 + linkDiameter))
            collisionPositions.pop_back();
            return true;
    }
    collisionPositions.pop_back();
    return false;
}