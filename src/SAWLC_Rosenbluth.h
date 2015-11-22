#ifndef H_SAWLC_ROSENBLUTH
#define H_SAWLC_ROSENBLUTH

#include "SAWLC.h"

using namespace std;

class SAWLC_Rosenbluth: public SAWLC
{
public:
	SAWLC_Rosenbluth(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_linkDiameter, double in_segConvFactor, 
                  double in_locPrecision, Eigen::Vector3d * in_initPoint);

	double getDefaultWeight();

	double getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where);

	double computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions,
                       int where);
};

class SACollector_Rosenbluth: public SACollector
{
public:
	vector<SAWLC_Rosenbluth*> myChains;
    SACollector_Rosenbluth(int in_numPaths,
             std::vector<double> & in_pathLength,
             std::vector<double> & in_linDensity,
             std::vector<double> & in_persisLength,
             std::vector<double> & in_linkDiameter,
             std::string in_nameDB,
             double in_segConvFactor = 1.0,
             double in_locPrecision = 0.0,
             bool in_fullSpecParam = false);

    SAWLC_Rosenbluth * getChainPointer(int i);
};

#endif