#include "Path.h"

extern std::default_random_engine randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

Path::Path(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_segConvFactor,
                  Eigen::Vector3d * in_initPoint)
{
    numPaths = in_numPaths;
    pathLength = in_pathLength;
    linDensity = in_linDensity;
    segConvFactor = in_segConvFactor;
    linkDiameter = 0.0;
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
    path.clear();                   //clears existing path
    path.push_back(initPoint);      //sets first point
    makePath(in_pathLength);        //creates the path
}

double Path::computeRg()
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

void Path::bumpPoints(double locPrecision)
{
    for (auto & oldPoint : path) {
        for (int i=0; i<3; i++) {
            oldPoint(i) += locPrecision * randNormalReal(randGenerator);
        }
    }
    return;
}

void Path::makePath(double in_pathLength) {
	cout << "Ran the wrong makePath!";
}



RgDict * Path::parSimChain()
{
    // work copies of variables
    vector<double> * Rg = new vector<double>(numPaths, 0.0);
    vector<double> * RgBump = new vector<double>(numPaths, 0.0);
    for (int i=0; i<numPaths; i++)
    {
        makeNewPath(pathLength.at(i));

        Rg->at(i) = computeRg();
        if (abs(locPrecision) > 0.001)
        {
            bumpPoints(locPrecision);
            RgBump->at(i) = computeRg();
        }
    }
    // we want to convert to user defined units - true
    RgDict * result = new RgDict(*Rg, *RgBump, 
        pathLength.at(0), linDensity,
        persisLength, linkDiameter, segConvFactor, true);
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




Collector::Collector(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 string in_nameDB,
                 double in_segConvFactor /*optional*/,
                 double in_locPrecision /*optional*/,
                 bool in_fullSpecParam /*optional*/)
{
	cout << "Calling base collector!" << endl;
	cout << "First linDensity element is: " << in_linDensity.at(0) << endl;
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

void Collector::startCollector()
{

    // Creates a vector of WLC entries, each vector entry is a
    // WLC instance pointer
    vector<RgDict*> data;
    // Compute the gyration radii for all the parameter pairs
    #pragma omp parallel for schedule(static,1) num_threads(4)
         // parallelize the next loop
        for (int i=0; i<linDensity.size(); i++) 
        {
            cout << "Simulation " << i << endl; cout.flush();
            data.push_back(getChainPointer(i)->parSimChain());
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
        dict->addToDBfileFull(fileDB); // adds each object
        cout << endl << "Density: " << dict->linDensity << ", Persistence length: "
             << dict->persisLength << std::endl << "Link diameter: "
             << dict->linkDiameter << std::endl;

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

}