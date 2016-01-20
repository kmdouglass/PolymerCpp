#include "Misc.h"

using namespace std;

unsigned seed;
std::minstd_rand randGenerator(1u);
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
    outVector.clear();
    for (int i=0; i<inVector.size(); i++)
    {
        if (multiplyBool) {
            outVector.push_back(inVector[i] * convFactor);
        }
        else {
            outVector.push_back(inVector[i] / convFactor);
        }
    }
    return;
}

double sum(vector<double> & vect)
{
    double result = 0.0;
    for (auto & val: vect)
        result += val;
    return result;
}
