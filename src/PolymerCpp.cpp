/* Classes for simulating random walk models for polymer
 * physics DNA/chromatin studies
 */

#include "PolymerCpp.h"


RgDict::RgDict(std::vector<double> * in_Rg, std::vector<double> * in_RgBump)
{
    Rg = *in_Rg;
    RgBump = *in_RgBump;
}


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

void Path::checkPath()
/* Checks that a path has the correct number of columns
 * NOT YET IMPLEMENTED
 */
{
    return;
}


WormlikeChain::WormlikeChain(double in_numSegments, int in_pLength,
                  Eigen::Vector3d * in_initPoint)
{
    numSegments = in_numSegments;
    pLength = in_pLength;
    initPoint = *in_initPoint;
    path.push_back(initPoint);
}

void WormlikeChain::makePath()
{
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

    // If numSegFrac is not zero, allocates memory also for leftover vector
    int numSeg = (abs(numSegFrac)>0.001) ? numSegInt : numSegInt-1;

    // Create the displacement distances in the tangent planes
    double *angDisp = new double[numSeg];
    double *tanPlaneDisp = new double[numSeg];
    double sigma = pow(2.0 / (double)pLength, 0.5);
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
        angDisp[numSeg-1] = pow(pLength*numSegFrac,0.5) *
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
    for (int i=0; i<numSeg; i++)
    {
        // Create a displacement in the plane tangent to currPoint
        dispVector = currPoint.cross(randVecs[i]);

        // Check if displacement and currPoint vectors are parallel
        while (dispVector.squaredNorm()<0.01)
        {
            Eigen::Vector3d * newRandVec = randPointSphere(1);
            dispVector = currPoint.cross(*newRandVec);
            delete newRandVec;
        }

        // Move the currPoint vector in the tangent plane
        dispVector = (dispVector / dispVector.squaredNorm()) *
                     tanPlaneDisp[i];

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
    delete angDisp, tanPlaneDisp;
    delete randVecs;
}

void WormlikeChain::makeNewPath()
{
    path.clear();               //clears existing path
    path.push_back(initPoint);  //sets first point
    makePath();                 //creates the path
    checkPath();                //checks if it is correct
}


Collector::Collector() { return; }

Collector::Collector(vector<double> & in_pathLength, 
              string in_nameDB,
              double in_segConvFactor /*optional*/)
{
    numPaths = in_pathLength.size();
    segConvFactor = in_segConvFactor;
    nameDB = in_nameDB;
    convSegments(pathLength, in_pathLength, segConvFactor, true);
}

void Collector::convSegments(vector<double> & outVector,
                      vector<double> & inVector,
                      double convFactor,
                      bool multiplyBool)
{
    for (int i=0; i<outVector.size(); i++)
    {
        if (multiplyBool) {
            outVector[i] = inVector[i] * convFactor;
        }
        else {
            outVector[i] = inVector[i] / convFactor;
        }
    }
    return;
}


WormlikeChainDict::WormlikeChainDict(WormlikeChain * in_chain,
                      double in_numSegments,
                      double in_locPrecision)
{
    chain = in_chain;
    numSegments = in_numSegments;
    locPrecision = in_locPrecision;
}

WormlikeChainDict::~WormlikeChainDict()
{
    delete chain;
}

WLCCollector::WLCCollector(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 string in_nameDB,
                 double in_segConvFactor /*optional*/,
                 double in_locPrecision /*optional*/,
                 bool in_fullSpecParam /*optional*/)
    {
        //copies variables to class
        pathLength = in_pathLength;
        linDensity = in_linDensity;
        numPaths = in_numPaths;
        persisLength = in_persisLength;
        nameDB = in_nameDB;
        segConvFactor = in_segConvFactor;
        locPrecision = in_locPrecision;
        fullSpecParam = in_fullSpecParam;

        // Convert from user-defined units to simulation units
        convSegments(pathLength, in_pathLength, segConvFactor, true);
        convSegments(linDensity, in_linDensity, segConvFactor, false);
        convSegments(persisLength, in_persisLength, segConvFactor, true);
        locPrecision = in_locPrecision * segConvFactor;

        // Starts the collector
        startCollector();
    }

void WLCCollector::startCollector()
{
    // Begin collecting wormlike chain conformation statistics.
    if (!fullSpecParam) 
    // Creates vectors with all permutations of linDensity and persisLength
    {
        vector<double> permLinDensity;
        vector<double> permPersisLength;
        for (int i=0; i<linDensity.size(); i++) {
            for (int j=0; j<persisLength.size(); j++) {
                permLinDensity.push_back(linDensity[i]);
                permPersisLength.push_back(persisLength[j]);
            }
        }
        linDensity = permLinDensity;
        persisLength = permPersisLength;
    }

    // Creates a vector of dictionaries, each vector entry is a
    // dictionary which has a WLC instance pointer and info about it
    vector<WormlikeChainDict*> myChains;
    for (int i=0; i<linDensity.size(); i++)
    {
        // get numSegments for different pathLengths for this linDensity
        vector<double> numSegments;
        for (int j=0; j<pathLength.size(); j++) {
            numSegments.push_back(pathLength[j]/linDensity[i]);
        }
        Eigen::Vector3d startDir(1.0,0.0,0.0);
        WormlikeChain * myChain = 
            new WormlikeChain(numSegments[0], persisLength[i],
                              &startDir);
        myChains.push_back(new WormlikeChainDict(
            myChain, numSegments[0], locPrecision));
    }
    vector<RgDict*> data;
    // Compute the gyration radii for all the parameter pairs
    for (int i=0; i<myChains.size(); i++) 
    {
        data.push_back(parSimChain(myChains[i]));
    }
    // SAVE THE CALCULATED RgData
}

RgDict * parSimChain(WormlikeChainDict * chainDict)
{
    // work copies of variables
    WormlikeChain * chain = chainDict->chain;
    float numSegments = chainDict->numSegments;
    float locPrecision = chainDict->locPrecision;

    //CHANGE to correct length later!
    int numPaths = 1;
    vector<double> * Rg = new vector<double>(numPaths, 0.0);
    vector<double> * RgBump = new vector<double>(numPaths, 0.0);
    Eigen::Vector3d randStartDir(0.0,0.0,0.0);
    for (int i=0; i<numPaths; i++)
    {
        for (int j=0; j<3; j++)
        {
            randStartDir(j) = randUniformReal(randGenerator) - 0.5;
        }
        chain->initPoint = randStartDir;
        chain->makeNewPath();

        Rg->at(i) = computeRg(chain, 3 /*dimensions*/);
        if (abs(locPrecision) > 0.001)
        {
            WormlikeChain * bumpedPath;
            bumpPoints(chain, locPrecision);
            RgBump->at(i) = computeRg(chain, 3 /*dimensions*/);
        }
    }
    RgDict * result = new RgDict(Rg, RgBump);
    return result;
}

double computeRg(WormlikeChain * chain, int dimensions /*optional*/)
{
    vector<Eigen::Vector3d> * path = &(chain->path);
    double Rg = 0.0;
    double secondMoments[3] = {0.0, 0.0, 0.0};
    double mean; double variance;
    for (int i=0; i<3; i++) {
        mean = 0.0; variance = 0.0;
        for (auto & point : *path ) { // iterate over all points in path
            mean += point(i);
        }
        mean = mean / path->size();

        for (auto & point : *path ) {
            variance += (mean-point(i))*(mean-point(i));
        }
        variance = variance / path->size();
        secondMoments[i] = variance;
    }

    if (dimensions == 3) {
        Rg = sqrt(secondMoments[0] + secondMoments[1] + secondMoments[2]);
    }
    else if (dimensions == 2) {
        Rg = sqrt(secondMoments[0] + secondMoments[1]);
    }
    return Rg;

}

void bumpPoints( WormlikeChain * chain, double locPrecision )
{
    vector<Eigen::Vector3d> * oldPath = &(chain->path);
    for (auto & oldPoint : *oldPath) {
        for (int i=0; i<3; i++) {
            oldPoint(i) += locPrecision * randNormalReal(randGenerator);
        }
    }
    return;
}



DateClass::DateClass() // constructor
{
    tick = std::time(&tick);
    timeinfo = localtime (&tick);
    char date[20], watch[20];
    strftime(date, 20, "%y-%m-%d", timeinfo);
    strftime(watch, 20, "%H-%M-%S", timeinfo);
    m_date.assign(date);
    m_time.assign(watch);
}
    /* redundant constructor, delete?
DateClass::DateClass(time_t tock)
{
    char date[20], watch[20];
    timeinfo = localtime (&tock);
    strftime(date, 20, "%y-%m-%d", timeinfo);
    strftime(watch, 20, "%H-%M-%S", timeinfo);
    m_date.assign(date);
    m_time.assign(watch);
}
    */
string DateClass::GetDate()
{
    return m_date;
}
string DateClass::GetTime()
{
    return m_time;
}

int main (int argc, char *argv[])
{
    DateClass DateTimeOnLaunch;  // Get time information for naming the database
    cout << "Today is " << DateTimeOnLaunch.GetDate() << ".\n";
    Eigen::Vector3d basePoint(8.,6.,5.);
    WormlikeChain WLC(3.,2.,&basePoint);
    cout << "Generating 10000 random sphere vectors.\n";
    WLC.points = WLC.randPointSphere(10000);
    for (int i=0; i<10000; i++)
    {
        printf("%f\t%f\t%f\n",WLC.points[i](0),WLC.points[i](1),
            WLC.points[i](2));
    }
    cout << "basePoint of WLC is: " << WLC.initPoint;
}
