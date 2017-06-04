#ifndef H_SAWLC_ROSENBLUTH
#define H_SAWLC_ROSENBLUTH

#include "SAWLC.h"

using namespace std;

class SAWLC_Rosenbluth: public SAWLC
{
public:
	SAWLC_Rosenbluth(double in_pathLength, 
                        double in_persisLength,
                        double in_linkDiameter,
                        Eigen::Vector3d * in_initPoint);
  // Constructor.

	double getDefaultWeight();
  /* There are many occasions, when there aren't any other points of chain
   * nearby, and the next step can only collide with the previous point.
   * All steps like these have the same weight, which we call defaultWeight,
   * so we can precompute it with higher precision and then just use the
   * cached value.
   */

	double getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where);
  /* Decides if you use defaultWeight or you compute the integral. */

	double computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions,
                       int where);
  /* Computes the step weight by carrying out a Monte Carlo integration
   * across the whole phase space. */
};

class SACollector_Rosenbluth: public SACollector
// Collector for Rosenbluth chains. Nothing surprising here.
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