#include "WLCCollector.h"

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