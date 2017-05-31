#include "SAWLC_Rosenbluth.h"

using namespace std;

extern std::minstd_rand randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

SAWLC_Rosenbluth::SAWLC_Rosenbluth(double in_pathLength, 
                                  double in_persisLength,
                                  double in_linkDiameter,
                                  Eigen::Vector3d * in_initPoint)
	: SAWLC(in_pathLength, in_persisLength, in_linkDiameter, in_initPoint)
{
	defaultWeight = getDefaultWeight();	
}

double SAWLC_Rosenbluth::getDefaultWeight()
{
    // Create simplest path configuration for calculation of this integral.
    vector<Eigen::Vector3d> basicPath(2, Eigen::Vector3d::Constant(0.0));
    basicPath[1](2)=1.0;
    Eigen::Vector3d nextPoint, currPoint;
    currPoint = basicPath[1];
    vector<int> collisions(1,0);

    double sigma = abs(persisLength)>0.00001 ? pow(2.0 / persisLength, 0.5) 
                                             : 0.0;

    Eigen::Vector3d dispVector;
    double projDistance;
    double angDisp, tanPlaneDisp;
    Eigen::Vector3d * randVecs = NULL;
    int passed = 0;
    randVecs = randPointSphere(30000);
    // See how many of these 30000 vectors are viable steps.
    for (int i=0; i<30000; i++)
    {
        if (sigma>0.00001) 
            angDisp = sigma * randNormalReal(randGenerator);
        else
            angDisp = 2*pi*randUniformReal(randGenerator);
        tanPlaneDisp = sin(angDisp);
        dispVector = randVecs[i].cross(currPoint);
        dispVector /= sqrt(dispVector.squaredNorm());
        dispVector = dispVector * tanPlaneDisp;
        projDistance = 1 - cos(angDisp);
        nextPoint = ((1 - projDistance) * currPoint) + dispVector;

        if (!checkCollision(&nextPoint, basicPath, collisions, 2))
            passed++;
    }
    delete randVecs;
    // Return the fraction of vectors that passed - this is the weight.
    return ((double)passed)/30000.0;
}

double SAWLC_Rosenbluth::computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions, int where)
{
    /* Similar as with defaultWeight, but this time we can have more possible
     * collision points, positions of which are stored in the collisionPositions
     * vector.
     */
    Eigen::Vector3d nextPoint, currPoint;
    double sigma = abs(persisLength)>0.00001 ? pow(2.0 / persisLength, 0.5) 
                                             : 0.0;
    Eigen::Vector3d dispVector;
    double projDistance;
    double angDisp, tanPlaneDisp;
    Eigen::Vector3d * randVecs = NULL;
    int passed = 0;
    randVecs = randPointSphere(5000);
    currPoint = cumulativePath[where-1]-cumulativePath[where-2];
    for (int i=0; i<5000; i++) // This function is done often, so less precision.
    {
        if (sigma>0.00001) 
            angDisp = sigma * randNormalReal(randGenerator);
        else
            angDisp = 2*pi*randUniformReal(randGenerator);
        tanPlaneDisp = sin(angDisp);
        dispVector = randVecs[i].cross(currPoint);
        dispVector /= sqrt(dispVector.squaredNorm());
        dispVector = dispVector * tanPlaneDisp;
        projDistance = 1 - cos(angDisp);
        nextPoint = ((1 - projDistance) * currPoint) + dispVector;
        if (!checkCollision(&nextPoint, cumulativePath, 
                            collisionPositions, where))
            passed++;
    }
    delete randVecs;
    return ((double)passed)/5000.0;  	
}

double SAWLC_Rosenbluth::getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where)
{
    /* If there is only one point available for collision, it's the previous one,
     * so the weight is the defaultWeight which we cleverly precomputed with
     * higher precision. */
    if (collisionPositions.size() == 1)
        return defaultWeight;
    /* If there are more possible collisions, we have to compute the integral.*/
    else
        return computeIntegral(cumulativePath, collisionPositions, where);
}



SACollector_Rosenbluth::SACollector_Rosenbluth(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 vector<double> & in_linkDiameter,
                 string in_nameDB,
                 double in_segConvFactor,
                 double in_locPrecision,
                 bool in_fullSpecParam)
: SACollector(in_numPaths, in_pathLength, in_linDensity, in_persisLength,
				in_linkDiameter,
                in_nameDB, in_segConvFactor, in_locPrecision, in_fullSpecParam)
{
    // Since it uses SACollector, it has to delete chains created by parent class
    // and then substitute its own chains.
    for (auto & chain: myChains)
    {
    	delete chain;
    }
    for (int i=0; i<linDensity.size(); i++)
    {
        Eigen::Vector3d startDir(1.0,0.0,0.0);
        SAWLC_Rosenbluth * myChain = 
            new SAWLC_Rosenbluth(pathLength[i], 
                                 persisLength[i], 
                                 linkDiameter[i],
                                 &startDir);

        myChains.push_back(myChain);
    }

}

SAWLC_Rosenbluth * SACollector_Rosenbluth::getChainPointer(int i)
{
	return myChains[i];
}
