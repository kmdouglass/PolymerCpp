#include "Path.h"

extern std::minstd_rand randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

Path::Path(double in_pathLength, 
          Eigen::Vector3d * in_initPoint)
{
    pathLength = in_pathLength;
    linkDiameter = 0.0; // default value, overriden in SAWLC
    chainWeight = 1.0; // default value, overriden in Rosenbluth
    initPoint = *in_initPoint;
    path.push_back(initPoint);
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


void Path::makeNewPath(double in_pathLength)
{
    // TODO: Account for fractional path lengths
    int numSegments = (int) in_pathLength;
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
    } while (path.size() < numSegments);
}

double Path::computeRg()
{
    double Rg = 0.0;
    double secondMoments[3] = {0.0, 0.0, 0.0};
    double mean; double variance;
    /* Radius of gyration is equal to sqrt of sums of variances
     * in each axis.
     */
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

void Path::bumpPoints(double locPrecision)
{
    for (auto & oldPoint : path) {
        for (int i=0; i<3/*EPFL!*/; i++) { // each point is bumped in each axis 
            oldPoint(i) += locPrecision * randNormalReal(randGenerator);
        }
    }
    return;
}

void Path::makePath(double in_pathLength) {
	cout << "Ran the wrong makePath!"; 
    // This obviously shouldn't happen, the function is virtual.
}



/* TODO: Needs to restructured like the rest of the project.
RgDict * Path::parSimChain()
{
    // work copies of variables
    vector<double> * Rg = new vector<double>(numPaths, 0.0);
    vector<double> * RgBump = new vector<double>(numPaths, 0.0);
    vector<double> * Wt = new vector<double>(numPaths, 0.0);
    double weightSum = 0.0;
    for (int i=0; i<numPaths; i++)
    // Create path, compute its gyration radius, weight, bump it
    // and add the values to respective vectors.
     
    {
        makeNewPath(pathLength.at(i));

        Rg->at(i) = computeRg();
        if (abs(locPrecision) > 0.001)
        {
            bumpPoints(locPrecision);
            RgBump->at(i) = computeRg();
        }
        Wt->at(i) = chainWeight;
        weightSum += Wt->at(i);
    }
    for (int i=0; i<numPaths; i++)
    // Normalize sum of all path weights to 1.
    {
        Wt->at(i) /= weightSum;
    }
    // we want to convert to user defined units - true
    RgDict * result = new RgDict(*Rg, *RgBump, *Wt,
        pathLength.at(0), linDensity,
        persisLength, linkDiameter, segConvFactor, true);
    return result;
}
*/

double theoreticalWLCRg(double Lp, double Lc)
{
    /* Computes the theoretical value for the chain's radius of gyration.
     * 
     * Parameters
     * ----------
     * Lp
     *     Persistence length of the chain.
     * Lc   
     *     Contour length of the chain.
     * 
     */
    double Rg2 = Lp * Lc / 3 - Lp*Lp + 2 * (Lp*Lp*Lp) / (Lc * Lc) *
                     (Lc - Lp * (1 - exp(-Lc/Lp)));
    return Rg2;
}




Collector::Collector(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 string in_nameDB,
                 double in_segConvFactor /*optional*/,
                 double in_locPrecision /*optional*/,
                 bool in_fullSpecParam /*optional*/)
{
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

    if (!fullSpecParam)
    // Create a grid of parameters.
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
}

Path * Collector::getChainPointer(int i) {}
// Virtual function.

void Collector::startCollector()
{
 
    /* TODO: Needs to restructured like the rest of the project.
    vector<RgDict*> data(linDensity.size()); // vector for storing data
    #pragma omp parallel for schedule(dynamic,1)
        // parallelize the next loop
        for (int i=0; i<linDensity.size(); i++) 
        {
            //cout << "Simulation " << i << endl; cout.flush();
            data.at(i)=getChainPointer(i)->parSimChain(); // most demanding part!
        }
    
    // Save data:

    // next line can be used to set maximum precision, but that should not be
    // necessary
    // fileDB << std::setprecision(std::numeric_limits<double>::digits10 + 1);
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
        dict->addToDBfileFull(fileDB); // adds each object
        /* SWITCH - just if we want to print data to check it is correct
        double sumRg = 0, sumRgBump = 0;
        for (int i=0; i<dict->Rg.size(); i++) {
            sumRg += dict->Rg.at(i)*dict->Wt.at(i);
            sumRgBump += dict->RgBump.at(i)*dict->Wt.at(i);
        }
        cout << "   Mean Rg:     " << sumRg << endl;
        cout << "   Mean RgBump: " << sumRgBump << " ";
        cout << "   Calculated:  " << theoreticalWLCRg(dict->linDensity,
                         dict->persisLength, dict->pathLength) << endl;
        //*/
        /*cout << dict->linkDiameter << " " 
             << dict->persisLength << " "
             << sumRg << std::endl; 
        //
    }
    fileDB.close();
*/
}
