#include "RgDict.h"

using namespace std;

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

std::vector<RgDict*> * readFromDBfile(vector<string> & dbNameList)
{
    // create new vector to hold dicts
    vector<RgDict*> * data;
    data = new vector<RgDict*>;

    for (auto & name: dbNameList)
    {
        std::ifstream fileDB;
        fileDB.open(name);
        if (!fileDB.is_open())
        {
            std::stringstream buffer;
            buffer << "Could not open file: " << name << std::endl;
            throw std::runtime_error(buffer.str());
        }
        // read vectorSize and segConvFactor
        int vectorSize; double segConvFactor;
        fileDB >> vectorSize;
        fileDB >> segConvFactor;


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
    }
    return data;
}

vector<RgDict*> * loadModel(vector<string> & dbNameList)
{
    vector<RgDict*> * data = new vector<RgDict*>;
    for (auto & name: dbNameList)
    {
        std::ifstream file;
        file.open(name);
        if (!file.is_open())
        {
            std::stringstream buffer;
            buffer << "Could not open file: " << name << std::endl;
            throw std::runtime_error(buffer.str());
        }
        int numObjects; double segConvFactor;
        file >> numObjects; file >> segConvFactor;
        
        for (int i=0; i<numObjects; i++)
        {

        }
    }
}
