

#include "WormlikeChain.h"

using namespace std;

extern std::default_random_engine randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

Eigen::Vector3d * Path::randPointSphere(int numPoints)
{
    Eigen::Vector3d * points = new Eigen::Vector3d[numPoints];
    for (int i=0; i<numPoints; i++) {
        
        // generate 2 random numbers and transform them into angles
        double phi = 2 * pi * randUniformReal(randGenerator);
        double theta = acos(2 * randUniformReal(randGenerator) - 1);
        
        // convert to cartesian coordinates and fill vector
        points[i](0) = cos(phi) * sin(theta);
        points[i](1) = sin(phi) * sin(theta);
        points[i](2) = cos(theta);
    }
    return points;
}

Eigen::Vector3d * Path::randPointSphere()
{
    Eigen::Vector3d * point = new Eigen::Vector3d;
	// generate 2 random numbers and transform them into angles
    double phi = 2 * pi * randUniformReal(randGenerator);
    double theta = acos(2 * randUniformReal(randGenerator) - 1);
        
    // convert to cartesian coordinates and fill vector
    point[0](0) = cos(phi) * sin(theta);
    point[0](1) = sin(phi) * sin(theta);
    point[0](2) = cos(theta);

    return point;
}


WormlikeChain::WormlikeChain(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_segConvFactor, double in_locPrecision, 
                  Eigen::Vector3d * in_initPoint)
{
    numPaths = in_numPaths;
    pathLength = in_pathLength;
    locPrecision = in_locPrecision;
    linDensity = in_linDensity;
    persisLength = in_persisLength;
    segConvFactor = in_segConvFactor;

    initPoint = *in_initPoint;
    path.push_back(initPoint);
}

void WormlikeChain::makePath(double in_pathLength)
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

    // Create the displacement distances in the tangent planes
    double *angDisp = new double[numSeg];
    double *tanPlaneDisp = new double[numSeg];
    double sigma = pow(2.0 / persisLength, 0.5);
    for (int i=0; i<numSeg; i++)
    {
        angDisp[i] = sigma * randNormalReal(randGenerator);
        tanPlaneDisp[i] = sin(angDisp[i]);
    }
    // Create random vectors uniformly sampled from the unit sphere
    Eigen::Vector3d * randVecs;
    randVecs = randPointSphere(numSeg); // don't forget to delete after

    // Final small displacement for non-integer numSegments
    if (abs(numSegFrac) > 0.001)
    {
        angDisp[numSeg-1] = pow(persisLength*numSegFrac,0.5) *
                       randNormalReal(randGenerator);
        tanPlaneDisp[numSeg-1] = numSegFrac * sin(angDisp[numSeg-1]);
        randVecs[numSeg-1] = numSegFrac * randVecs[numSeg-1];
    }

    // Primary iterative loop for creating the chain
    Eigen::Vector3d currPoint = initPoint;
    Eigen::Vector3d * workingPath = new Eigen::Vector3d[numSeg];
    Eigen::Vector3d dispVector, nextPoint;
    double projDistance;
    workingPath[0] = currPoint;

    for (int i=1; i<numSeg; i++)
    {
        // Create a displacement in the plane tangent to currPoint
        dispVector = currPoint.cross(randVecs[i]);
        dispVector /= sqrt(dispVector.squaredNorm());
        // Check if displacement and currPoint vectors are parallel
        int j=0;
        while (dispVector.squaredNorm()<0.01)
        {
            currPoint /= sqrt(currPoint.squaredNorm());
            cout << "while   " << j++; cout.flush();
            Eigen::Vector3d * newRandVec = randPointSphere(1);
            dispVector = currPoint.cross(*newRandVec);
            cout << *newRandVec << endl;
            cout << currPoint << endl;
            cout << dispVector << endl;
            cout << "Norm: " << dispVector.squaredNorm() << endl;
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
    for (int i=1; i<numSeg; i++)
    {
        path.push_back(workingPath[i] + path[i-1]);
    }
    delete workingPath;
    delete angDisp; delete tanPlaneDisp;
    delete randVecs;
}

void WormlikeChain::makeNewPath(double in_pathLength)
{
    path.clear();                   //clears existing path
    path.push_back(initPoint);      //sets first point
    makePath(in_pathLength);        //creates the path
}

double WormlikeChain::computeRg()
{
    double Rg = 0.0;
    double secondMoments[3] = {0.0, 0.0, 0.0};
    double mean; double variance;
    for (int i=0; i<3; i++) {
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


    Rg = sqrt(secondMoments[0] + secondMoments[1] + secondMoments[2]);
    return Rg;

}

void WormlikeChain::bumpPoints(double locPrecision)
{
    for (auto & oldPoint : path) {
        for (int i=0; i<3; i++) {
            oldPoint(i) += locPrecision * randNormalReal(randGenerator);
        }
    }
    return;
}

RgDict * parSimChain(WormlikeChain * chain)
{
    // work copies of variables
    vector<double> * Rg = new vector<double>(chain->numPaths, 0.0);
    vector<double> * RgBump = new vector<double>(chain->numPaths, 0.0);
    for (int i=0; i<chain->numPaths; i++)
    {
        for (int j=0; j<3; j++)
        {
            chain->initPoint(j) = randUniformReal(randGenerator) - 0.5;
            // normalize, otherwise can be too small
            chain->initPoint /= sqrt(chain->initPoint.squaredNorm());
        }
        chain->makeNewPath(chain->pathLength.at(i));

        Rg->at(i) = chain->computeRg();
        if (abs(chain->locPrecision) > 0.001)
        {
            chain->bumpPoints(chain->locPrecision);
            RgBump->at(i) = chain->computeRg();
        }
    }
    // we want to convert to user defined units - true
    RgDict * result = new RgDict(Rg, RgBump, 
        chain->pathLength.at(0), chain->linDensity,
        chain->persisLength, chain->segConvFactor, true); 
    return result;
}

double theoreticalWLCRg(double c, double Lp, double N)
{
    double Rg2 = (Lp * N / c) / 3 - 
                    Lp*Lp + 
                    2 * (Lp*Lp*Lp) / ((N / c)*(N / c)) *
                  ((N / c) - Lp * (1 - exp(-(N/c)/Lp)));
    return sqrt(Rg2);
}


