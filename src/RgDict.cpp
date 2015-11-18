#include "RgDict.h"

using namespace std;

RgDict::RgDict(std::vector<double> & in_Rg, std::vector<double> & in_RgBump,
    double in_pathLength, double in_linDensity, double in_persisLength, 
    double in_linkDiameter, double in_segConvFactor, bool convert)
{
    Rg = in_Rg;
    RgBump = in_RgBump;
    pathLength = in_pathLength;
    if (convert) {
        convSegments(Rg, in_Rg, in_segConvFactor, false);
        convSegments(RgBump, in_RgBump, in_segConvFactor, false);
        linDensity = in_linDensity * in_segConvFactor;
        persisLength = in_persisLength / in_segConvFactor;
        linkDiameter = in_linkDiameter / in_segConvFactor;
    }
    else {
        linDensity = in_linDensity;
        persisLength = in_persisLength;
        linkDiameter = in_linkDiameter;
    }
    segConvFactor = in_segConvFactor;
}

void RgDict::addToDBfileFull(std::ofstream & fileDB)
{
    fileDB << endl << pathLength;
    fileDB << " " << linDensity;
    fileDB << " " << persisLength;
    fileDB << " " << linkDiameter;
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
    fileDB << endl;
}

void RgDict::addToDBfileShort(std::ofstream & fileDB)
{
    fileDB << endl << pathLength;
    fileDB << " " << linDensity;
    fileDB << " " << persisLength;
    fileDB << " " << linkDiameter;
    fileDB << " " << Rg.size() << endl;
    double mean = 0.0;
    for (int i=0; i<Rg.size(); i++)
    {
        mean += Rg[i];
        
    }
    mean /= (double)Rg.size();
    fileDB << " " << mean;

    double meanBump = 0.0;
    for (int i=0; i<Rg.size(); i++)
    {
        meanBump += RgBump[i];
    }
    meanBump /= (double)Rg.size();
    fileDB << " " << meanBump << std::endl;
}
