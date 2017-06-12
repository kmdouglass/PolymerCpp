#include "SAWLC.h"

using namespace std;

extern std::minstd_rand randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

SAWLC::SAWLC(double in_pathLength, 
              double in_persisLength,
              double in_linkDiameter, 
              Eigen::Vector3d * in_initPoint)
    : Path(in_pathLength, in_initPoint)
{
    persisLength = in_persisLength;
    linkDiameter = in_linkDiameter;
    defaultWeight = 1.0;

    if (1.0 < linkDiameter)
    {
        std::stringstream buffer;
        buffer << "The link diameter is larger than the distance between links. "
            << std::endl << "Link diameter: " << linkDiameter << std::endl
            << "Distance between links: " << 1.0 << std::endl;
        throw std::out_of_range(buffer.str());
    }
}

void SAWLC::makePath(double in_pathLength)
{   
    /* this chunk of code is the same in SAWLC and WLC
     * but I cant put it to Path:: because it would go out of scope
     */
    
    // Add one because two vertices = one segment
    int numVerts = (int) in_pathLength + 1;
    // check if numVerts is in valid range
    if (numVerts <3)
    {
        std::stringstream buffer;
        buffer << "The number of segments must be greater than 3, but a "
            << "value of " << numVerts << "was supplied. "
            << "Why would you want such a short chain anyways?" << std::endl;
        throw std::out_of_range(buffer.str());
    }
    Eigen::Vector3d currPoint = initPoint;
    Eigen::Vector3d dispVector, nextPoint;
    double projDistance;

    // Primary iterative loop for creating the chain
    vector<Eigen::Vector3d> workingPath(numVerts, Eigen::Vector3d::Constant(0.0));
    vector<Eigen::Vector3d> cumulativePath(numVerts, Eigen::Vector3d::Constant(0.0));
    workingPath[1] = initPoint;
    cumulativePath[1] = initPoint;
    vector<double> stepWeights(numVerts,0.0);
    stepWeights[0] = 1.0; stepWeights[1] = defaultWeight;
    vector<int> collisionPositions;
    double angDisp, tanPlaneDisp;
    Eigen::Vector3d * randVec = NULL;
    for (int i=2; i<numVerts; i++)
    {
        // Scan the chain for possible collisions and save them to colPos vector
        getPossibleCollisions(cumulativePath, collisionPositions, i);

        // Calculate step weight (in SAWLC it is always 1)
        stepWeights[i] = getNextStepWeight(cumulativePath, collisionPositions, i);
        int j = 0;
        while (true)
        // keep generating vectors until you get one which is ok
        {
            j++;
            
            angDisp = pow(-2.0 / persisLength * log(1 - randUniformReal(randGenerator)), 0.5);
            tanPlaneDisp = sin(angDisp);
            
            randVec = randPointSphere();
            dispVector = randVec->cross(currPoint);
            dispVector /= sqrt(dispVector.squaredNorm());
            dispVector = dispVector * tanPlaneDisp;
            projDistance = 1 - cos(angDisp);
            nextPoint = ((1 - projDistance) * currPoint) + dispVector;
            if (!checkCollision(&nextPoint, cumulativePath, collisionPositions, i))
            // step is OK, add it to path and move on
            {
                workingPath[i] = nextPoint;;
                cumulativePath[i] = nextPoint + cumulativePath[i-1];
                currPoint = nextPoint;
                delete randVec;
                break;
            }
            // Otherwise, delete newly allocated vectors and try again.
            delete randVec;
            // If 10000 attempts fail, throws exception which resets the chain.
            if (j%10000 == 0)
            {
                //cout << "DEAD END" << std::endl; cout.flush();
                return;
            }
        }
    	// move to next step
    }

    // copy final vector to path
    path.clear();
    for (int i=0; i<numVerts; i++)
    {
        path.push_back(cumulativePath[i]);
    }

    // calculate and save weight
    double wt = 1.0;
    for (auto & entry: stepWeights)
        wt *= entry;
    chainWeight=wt;
}

double SAWLC::getDefaultWeight()
{
    return 1.0;
}

double SAWLC::computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions, int where)
{
    return 1.0;
}

double SAWLC::getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where)
{
    return 1.0;
}

void SAWLC::getPossibleCollisions(
                    vector<Eigen::Vector3d> & cumulativePath,
                    vector<int> & collisionPositions,
                    int where)
{
    collisionPositions.clear();
    Eigen::Vector3d distanceVector; double distanceSquared;
    Eigen::Vector3d endPoint = cumulativePath.at(where-1);
    // we use square of distance, because we dont want to compute sqrt()
    // in each iteration, and distance is always positive
    // in each step this is the furthest that the collision sphere can reach
    double limitSquared = (1.0 + linkDiameter)*(1.0 + linkDiameter);
    for (int i=0; i<where-2; i++)
    // ignore last 2 points because they are always a potential collision
    {
        // get distance of this point from end of chain
        distanceVector = endPoint - cumulativePath[i];
        distanceSquared = distanceVector.squaredNorm();
        if (distanceSquared < limitSquared)
        // point is close enough that it might cause a collision
        {
            collisionPositions.push_back(i);
            /*cout << "Possible collision found at link no. "
                << where << ", with link no. "
                << i << std::endl; cout.flush();*/
        }
    }
    collisionPositions.push_back(where-2);
    // the second-to-last link is always available for collision
}

bool SAWLC::checkCollision(Eigen::Vector3d * nextPoint,
                           vector<Eigen::Vector3d> & cumulativePath,
                           vector<int> & collisionPositions,
                           int where)
{
    Eigen::Vector3d endPoint = *nextPoint + cumulativePath[where-1];
    Eigen::Vector3d distanceVector; double distanceSquared;
    double limitSquared = linkDiameter*linkDiameter;

    for (auto & position : collisionPositions)
    /* Scans every point in collisionPositions (which was acquired beforehand)
     * and checks if the step collides with it */
    {
        distanceVector = endPoint - cumulativePath[position];
        distanceSquared = distanceVector.squaredNorm();
        if (distanceSquared < limitSquared)
        {
            /*cout << "Collision averted at link no. "
                << where-1 << std::endl
                << "Distance: " << sqrt(distanceSquared) << std::endl
                << "linkDiameter: " << linkDiameter << std::endl; cout.flush();
            */
            return true;
        }
    }
    return false;
}


SACollector::SACollector(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 vector<double> & in_linkDiameter,
                 string in_nameDB,
                 double in_segConvFactor,
                 double in_locPrecision,
                 bool in_fullSpecParam)
    : Collector(in_numPaths, in_pathLength, in_linDensity, in_persisLength,
                in_nameDB, in_segConvFactor, in_locPrecision, in_fullSpecParam)
{
    linkDiameter = in_linkDiameter;
    convSegments(linkDiameter, in_linkDiameter, segConvFactor, true);
    //copies variables to class
    if (!fullSpecParam)
    {
        vector<double> permLinkDiameter;
        vector<double> permLinDensity;
        vector<double> permPersisLength;
        for (int i=0; i<linDensity.size(); i++) {
            for (int j=0; j<linkDiameter.size(); j++) {
                permLinDensity.push_back(linDensity[i]);
                permPersisLength.push_back(persisLength[i]);
                permLinkDiameter.push_back(linkDiameter[j]);
            }
        }
        linDensity = permLinDensity;
        persisLength = permPersisLength;
        linkDiameter = permLinkDiameter;
    }

    // Initializes chains, which are then simulated by parSimChain
    for (int i=0; i<linDensity.size(); i++)
    {
        Eigen::Vector3d startDir(1.0,0.0,0.0);
        SAWLC * myChain = 
            new SAWLC(pathLength[i],
                      persisLength[i],
                      linkDiameter[i], 
                      &startDir);

        myChains.push_back(myChain);
    }

}

SAWLC * SACollector::getChainPointer(int i)
{
    return myChains[i];
}

SACollector::~SACollector()
{
    for (auto & chain: myChains)
    {
        delete chain;
    }
}
