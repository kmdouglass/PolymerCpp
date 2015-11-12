#include "Misc.h"

using namespace std;

unsigned seed;
std::default_random_engine randGenerator(1u);
std::uniform_real_distribution<double> randUniformReal(0.0,1.0);
std::normal_distribution<double> randNormalReal(0.0,1.0);

void seedRandom()
{
    seed = std::chrono::system_clock::now().time_since_epoch().count();
    randGenerator.seed(seed);
}

void printRandom()
{
	cout << "Random number incoming: " << randNormalReal(randGenerator)
		 << std::endl;
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