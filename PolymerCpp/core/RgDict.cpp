#include "RgDict.h"

using namespace std;

RgDict::RgDict(std::vector<double> & in_Rg, std::vector<double> & in_RgBump,
    std::vector<double> & in_Wt,
    double in_pathLength, double in_linDensity, double in_persisLength, 
    double in_linkDiameter, double in_segConvFactor, bool convert)
{
    Rg = in_Rg;
    RgBump = in_RgBump;
    Wt = in_Wt;
    pathLength = in_pathLength;
    if (convert) {
        /* If we want to convert segments back to user units,
         * perform the necessary conversions. */ 
        convSegments(Rg, in_Rg, in_segConvFactor, false);
        convSegments(RgBump, in_RgBump, in_segConvFactor, false);
        linDensity = in_linDensity * in_segConvFactor;
        persisLength = in_persisLength / in_segConvFactor;
        linkDiameter = in_linkDiameter / in_segConvFactor;
    }
    else {
        /* Otherwise just copy variables. */
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
    /* SWITCH if you want weights in output
    for (int i=0; i<Rg.size(); i++)
    {
        fileDB << " " << Wt[i];
    } //*/
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
    // Shortened version of output - not used at the moment.
    fileDB << endl << pathLength;
    fileDB << " " << linDensity;
    fileDB << " " << persisLength;
    fileDB << " " << linkDiameter;
    fileDB << " " << Rg.size() << endl;
    double mean = 0.0;
    for (int i=0; i<Rg.size(); i++)
    {
        mean += Rg[i]*Wt[i];
        
    }
    mean /= (double)Rg.size();
    fileDB << " " << mean;

    double meanBump = 0.0;
    for (int i=0; i<Rg.size(); i++)
    {
        meanBump += RgBump[i]*Wt[i];
    }
    meanBump /= (double)Rg.size();
    fileDB << " " << meanBump << std::endl;
}

double RgDict::getVariance(bool bumped)
{
    // Expecting sum of all weights to be 1
    assert((abs(sum(Wt))<1.02 && abs(sum(Wt))>0.98));
    
    double variance = 0.0, meanOfSquare = 0.0, squareOfMean = 0.0;
    vector<double> * RgPtr;
    // Decide if we want variance of bumped chain or original one.
    RgPtr = bumped ? &RgBump : &Rg;

    double value;
    for (int i=0; i<RgPtr->size(); i++)
    {
        value = RgPtr->at(i);
        meanOfSquare += Wt.at(i) * value*value;
        squareOfMean += Wt.at(i) * value;
    }
    squareOfMean = squareOfMean * squareOfMean;

    variance = meanOfSquare - squareOfMean;
    return variance;
}