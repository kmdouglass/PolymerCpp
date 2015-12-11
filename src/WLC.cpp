#include "WLC.h"

using namespace std;

extern std::default_random_engine randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;


WLC::WLC(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_segConvFactor, double in_locPrecision, 
                  Eigen::Vector3d * in_initPoint) 
    : Path(in_numPaths, in_pathLength, in_linDensity,
            in_segConvFactor, in_initPoint)
{
    locPrecision = in_locPrecision;
    persisLength = in_persisLength;
}

void WLC::makePath(double in_pathLength)
{
    /* this chunk of code is the same in SAWLC and WLC
     * but I cant put it to Path:: because it would go out of scope
     */
    int numSegments = (int) (in_pathLength / linDensity);
    // check if numSegments is in valid range
    if (numSegments < 1)
    {
        std::stringstream buffer;
        buffer << "The number of segments must be greater than 3, but a "
            << "value of " << numSegments << "was supplied. "
            << "Why would you want such a short chain anyways?" << std::endl;
        throw std::out_of_range(buffer.str());
    }

    // split numSegments into integer and leftover
    // no leftover for the time being
    int numSeg;
    numSeg = (int)numSegments;

    double sigma = abs(persisLength)>0.00001 ? pow(2.0 / persisLength, 0.5) 
                                             : 0.0;
    Eigen::Vector3d currPoint = initPoint;
    Eigen::Vector3d dispVector, nextPoint;
    double projDistance;

    // Create the displacement distances in the tangent planes
    double *angDisp = new double[numSegments];
    double *tanPlaneDisp = new double[numSegments];

    //since we dont have to check for collisions,
    //we can pre-generate all required random numbers
    for (int i=0; i<numSegments; i++)
    {
        if (sigma>0.001) 
            angDisp[i] = sigma * randNormalReal(randGenerator);
        else
            angDisp[i] = 2*pi*randUniformReal(randGenerator);
        tanPlaneDisp[i] = sin(angDisp[i]);
    }
    // Create random vectors uniformly sampled from the unit sphere
    Eigen::Vector3d * randVecs;
    randVecs = randPointSphere(numSegments); // don't forget to delete after

    // Primary iterative loop for creating the chain
    Eigen::Vector3d * workingPath = new Eigen::Vector3d[numSegments];
    for (int i=1; i<numSegments; i++)
    {
        // Create a displacement in the plane tangent to currPoint
        dispVector = currPoint.cross(randVecs[i]);
        dispVector /= sqrt(dispVector.squaredNorm());

        // Check if displacement and currPoint vectors are parallel
        // This should not happen (you would have to be extremely)
        // lucky to get exactly parallel vectors with double precision
        int j=0;
        while (dispVector.squaredNorm()<0.01)
        {
            currPoint /= sqrt(currPoint.squaredNorm());
            Eigen::Vector3d * newRandVec = randPointSphere(1);
            dispVector = currPoint.cross(*newRandVec);
            delete newRandVec;
        }

        // Move the currPoint vector in the tangent plane
        dispVector = dispVector * tanPlaneDisp[i];

        // Back Project new point onto sphere
        projDistance = 1 - cos(angDisp[i]);
        nextPoint = ((1 - projDistance) * currPoint) + dispVector;

        // Append nextPoint to array of vectors on the path
        workingPath[i] = nextPoint;
        currPoint = nextPoint;
    }

    // Add up the vectors in path to create the polymer
    path[0] = workingPath[0];
    for (int i=1; i<numSegments; i++)
    {
        path.push_back(workingPath[i] + path[i-1]);
    }
    delete workingPath;
    delete angDisp; delete tanPlaneDisp;
    delete randVecs;
}

double WLC::getTheoreticalRg()
{
    return theoreticalWLCRg(linDensity, persisLength, pathLength[0]);
}


WLCCollector::WLCCollector(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 string in_nameDB,
                 double in_segConvFactor /*optional*/,
                 double in_locPrecision /*optional*/,
                 bool in_fullSpecParam /*optional*/)
    : Collector(in_numPaths, in_pathLength, in_linDensity, in_persisLength,
                in_nameDB, in_segConvFactor, in_locPrecision, in_fullSpecParam)
{
    for (int i=0; i<linDensity.size(); i++)
    {
        Eigen::Vector3d startDir(1.0,0.0,0.0);
        WLC * myChain = 
            new WLC(numPaths, pathLength, linDensity[i], 
                    persisLength[i], segConvFactor, locPrecision,
                    &startDir);
        myChains.push_back(myChain);
    }
}

WLC * WLCCollector::getChainPointer(int i)
{
    return myChains[i];
}

WLCCollector::~WLCCollector()
{
    for (auto & chain: myChains)
    {
        delete chain;
    }
}