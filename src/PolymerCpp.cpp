/* Classes for simulating random walk models for polymer
 * physics DNA/chromatin studies
 */


#include "PolymerCpp.h"


RgDict::RgDict(std::vector<double> * in_Rg, std::vector<double> * in_RgBump,
    double in_pathLength, double in_linDensity, double in_persisLength, 
    double in_segConvFactor, bool convert)
{
    Rg = *in_Rg;
    RgBump = *in_RgBump;
    pathLength = in_pathLength;
    if (convert) {
        convSegments(Rg, Rg, in_segConvFactor, false);
        convSegments(RgBump, RgBump, in_segConvFactor, false);
        linDensity = in_linDensity * in_segConvFactor;
        persisLength = in_persisLength / in_segConvFactor;
    }
    else {
        linDensity = in_linDensity;
        persisLength = in_persisLength;
    }
    segConvFactor = in_segConvFactor;
}

void RgDict::addToDBfile(std::ofstream & fileDB)
{
    fileDB << endl << pathLength;
    fileDB << " " << linDensity;
    fileDB << " " << persisLength;
    fileDB << " " << Rg.size() << endl;
    for (int i=0; i<Rg.size(); i++)
    {
        fileDB << " " << Rg[i];
    }
    fileDB << endl;
    for (int i=0; i<Rg.size(); i++)
    {
        fileDB << " " << RgBump[i];
    }
}

std::vector<RgDict*> * readFromDBfile(std::ifstream & fileDB)
{
    // read vectorSize and segConvFactor
    int vectorSize; double segConvFactor;
    fileDB >> vectorSize;
    fileDB >> segConvFactor;

    // create new vector to hold dicts
    vector<RgDict*> * data;
    data = new vector<RgDict*>;
    // read each vector in file
    for (int i=0; i<vectorSize; i++)
    {
        // load all values
        double linDensity; double persisLength; double pathLength;
        fileDB >> pathLength;
        fileDB >> linDensity;
        fileDB >> persisLength;
        int dictSize;
        fileDB >> dictSize;
        vector<double> Rg; vector<double> RgBump;
        double val;
        for (int j=0; j<dictSize; j++)
        {
            fileDB >> val;
            Rg.push_back(val);
        }
        for (int j=0; j<dictSize; j++)
        {
            fileDB >> val;
            RgBump.push_back(val);
        }
        // already in user-defined units - false
        RgDict * dict = new RgDict(&Rg, &RgBump, pathLength, 
            linDensity, persisLength, segConvFactor, false);
        data->push_back(dict);
    }
    return data;
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
    checkPath();                    //checks if it is correct
}

/*
WormlikeChainDict::WormlikeChainDict(WormlikeChain * in_chain,
                      vector<double> * in_numSegments,
                      double in_locPrecision,
                      double in_linDensity,
                      double in_persisLength,
                      double in_segConvFactor)
{
    chain = in_chain;
    numSegments = *in_numSegments;
    locPrecision = in_locPrecision;
    linDensity = in_linDensity;
    persisLength = in_persisLength;
    segConvFactor = in_segConvFactor;
}

WormlikeChainDict::~WormlikeChainDict()
{
    delete chain;
}
*/

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
        convSegments(linDensity, in_linDensity, segConvFactor, false);
        convSegments(persisLength, in_persisLength, segConvFactor, true);
        locPrecision = in_locPrecision * segConvFactor;


        // Starts the collector
        startCollector();
    }

void convSegments(vector<double> & outVector,
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

void WLCCollector::startCollector()
{
    cout << "Starting collector!" << endl;
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
    vector<WormlikeChain*> myChains;
    for (int i=0; i<linDensity.size(); i++)
    {
        Eigen::Vector3d startDir(1.0,0.0,0.0);
        WormlikeChain * myChain = 
            new WormlikeChain(numPaths, pathLength, linDensity[i], 
                    persisLength[i], segConvFactor, locPrecision,
                    &startDir);
        myChains.push_back(myChain);
    }
    vector<RgDict*> data;
    // Compute the gyration radii for all the parameter pairs
    #pragma omp parallel for // parallelize the next loop
    for (int i=0; i<myChains.size(); i++) 
    {
        cout << "Simulation " << i << endl; cout.flush();
        data.push_back(parSimChain(myChains[i]));
    }
    
    
    // SAVE THE CALCULATED RgData
    // Opens or creates the database file for adding simulation results

    // sets maximum output precision
    fileDB << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    fileDB.open(nameDB + ".txt", ios::out | ios::trunc);
    if (!fileDB.is_open())
    {
        std::stringstream buffer;
        buffer << "Could not open file: " << nameDB << ".txt" << std::endl;
        throw std::runtime_error(buffer.str());
    }


    // adds info about number of RgData objects to be saved
    fileDB << data.size();
    fileDB << " " << segConvFactor; // info about conversion factor

    for (auto & dict: data)
    {
        dict->addToDBfile(fileDB); // adds each object
        cout << "Density: " << dict->linDensity << ", Persistence length: "
             << dict->persisLength << std::endl;

        double sumRg = 0, sumRgBump = 0;
        for (int i=0; i<dict->Rg.size(); i++) {
            sumRg += dict->Rg.at(i);
            sumRgBump += dict->RgBump.at(i);
        }
        sumRg /= dict->Rg.size(); sumRgBump /= dict->RgBump.size();
        cout << "   Mean Rg:     " << sumRg << endl;
        cout << "   Mean RgBump: " << sumRgBump << endl;
        cout << "   Calculated:  " << theoreticalWLCRg(dict->linDensity,
                         dict->persisLength, dict->pathLength) << endl;
    }
    fileDB.close();


    // collect garbage
    for (auto & chain: myChains) delete chain;

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

        Rg->at(i) = chain->computeRg(3 /*dimensions*/);
        if (abs(chain->locPrecision) > 0.001)
        {
            chain->bumpPoints(chain->locPrecision);
            RgBump->at(i) = chain->computeRg(3 /*dimensions*/);
        }
    }
    // we want to convert to user defined units - true
    RgDict * result = new RgDict(Rg, RgBump, 
        chain->pathLength.at(0), chain->linDensity,
        chain->persisLength, chain->segConvFactor, true); 
    return result;
}

double WormlikeChain::computeRg(int dimensions /*optional*/)
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

    if (dimensions == 3) {
        Rg = sqrt(secondMoments[0] + secondMoments[1] + secondMoments[2]);
    }
    else if (dimensions == 2) {
        Rg = sqrt(secondMoments[0] + secondMoments[1]);
    }
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

double theoreticalWLCRg(double c, double Lp, double N)
{
    double Rg2 = (Lp * N / c) / 3 - Lp*Lp + 2 * (Lp*Lp*Lp) / ((N / c)*(N / c)) *
                  ((N / c) * Lp * (1 - exp(-(N/c)/Lp)));
    return sqrt(Rg2);
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

void seedRandom()
{
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    randGenerator.seed(seed);
}

int main (int argc, char *argv[])
{
    seedRandom();
    // Test case 1: Print time, test sphere sampling
    /*
    DateClass DateTimeOnLaunch;  // Get time information for naming the database
    cout << "Today is " << DateTimeOnLaunch.GetDate() << ".\n";

    Eigen::Vector3d startPoint(0.0,1.0,0.0);
    WormlikeChain WLC(1.,1.,&startPoint);
    cout << "Generating 10000 random sphere vectors.\n";
    WLC.points = WLC.randPointSphere(10000);

    double meanX = 0.0, meanY = 0.0, meanZ = 0.0, meanSquare = 0.0;
    for (int i=0; i<10000; i++)
    {
        sumX += WLC.points[i](0);
        sumY += WLC.points[i](1);
        sumZ += WLC.points[i](2);
        sumSquare += WLC.points[i](0)*WLC.points[i](0) +
                     WLC.points[i](1)*WLC.points[i](1) +
                     WLC.points[i](2)*WLC.points[i](2);
    }
    cout << "Sum of X,Y,Z: " << sumX << ", " << sumY
         << ", " << sumZ << std::endl;
    cout << "Sum of squares of vecors: " << sumSquare << std::endl;
    */

    //Test case 2: test for vector normalization
    /*
    Eigen::Vector3d point(6.3,0.7,-5.6);
    cout << "Vector: " << point << endl;
    cout << "Squared norm: " << point.squaredNorm() << endl;
    point /= sqrt(point.squaredNorm());
    cout << "Normalized vector: " << point << endl;
    cout << "Squared norm: " << point.squaredNorm() << endl;
    */

    // Test case 3: create a single random walk
    /*
    Eigen::Vector3d startPoint(0.0,1.0,0.0);
    WormlikeChain WLC(10.05, 25, &startPoint);
    WLC.makeNewPath();
    for (auto & point: WLC.path)
    {
        cout << point(0) << ", " << point(1) << ", " << point(2) << endl;
    }
    */

    // Test case 4: Create a WLCCollector,
    //              see if computed Rg matches theory
    
    
    int numPaths = 1000;

    vector<double> pathLength;
    for (int i=0; i<numPaths; i++)
        pathLength.push_back(25000);

    double linDensities[] = {10.0};
    vector<double> linDensity (linDensities, linDensities 
                                + sizeof(linDensities) / sizeof(double) );

    double persisLengths[] = {100.0};
    vector<double> persisLength (persisLengths, persisLengths 
                                + sizeof(persisLengths) / sizeof(double) );
    double segConvFactor = 2.5;
    WLCCollector myCollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             "data",
                             segConvFactor,
                             2.12,
                             true);
    
}