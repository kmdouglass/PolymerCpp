/* Classes for simulating random walk models for polymer
 * physics DNA/chromatin studies
 */


#include "Main.h"

extern std::default_random_engine randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;



int main (int argc, char *argv[])
{
    seedRandom();
    // Test case 1: Print time, test sphere sampling
    /*
    DateClass DateTimeOnLaunch;  // Get time information for naming the database
    cout << "Today is " << DateTimeOnLaunch.GetDate() << ".\n";

    Eigen::Vector3d startPoint(0.0,1.0,0.0);
    WormlikeChain WLC(1.,1.,&startPoint);
    cout << "Generating 10000 random sphere vectors.\n";
    WLC.points = WLC.randPointSphere(10000);

    double meanX = 0.0, meanY = 0.0, meanZ = 0.0, meanSquare = 0.0;
    for (int i=0; i<10000; i++)
    {
        sumX += WLC.points[i](0);
        sumY += WLC.points[i](1);
        sumZ += WLC.points[i](2);
        sumSquare += WLC.points[i](0)*WLC.points[i](0) +
                     WLC.points[i](1)*WLC.points[i](1) +
                     WLC.points[i](2)*WLC.points[i](2);
    }
    cout << "Sum of X,Y,Z: " << sumX << ", " << sumY
         << ", " << sumZ << std::endl;
    cout << "Sum of squares of vecors: " << sumSquare << std::endl;
    */

    //Test case 2: test for vector normalization
    /*
    Eigen::Vector3d point(6.3,0.7,-5.6);
    cout << "Vector: " << point << endl;
    cout << "Squared norm: " << point.squaredNorm() << endl;
    point /= sqrt(point.squaredNorm());
    cout << "Normalized vector: " << point << endl;
    cout << "Squared norm: " << point.squaredNorm() << endl;
    */

    // Test case 3: create a single random walk
    /*
    Eigen::Vector3d startPoint(0.0,1.0,0.0);
    WormlikeChain WLC(10.05, 25, &startPoint);
    WLC.makeNewPath();
    for (auto & point: WLC.path)
    {
        cout << point(0) << ", " << point(1) << ", " << point(2) << endl;
    }
    */

    // Test case 4: Create a WLCCollector,
    //              see if computed Rg matches theory
    
    int numPaths; double input;
    cout << "Enter numPaths: ";
    cin >> numPaths;
    vector<double> pathLength;
    double in_pathLength;
    cout << "Enter base pathLength: ";
    cin >> in_pathLength;
    for (int i=0; i<numPaths; i++)
        pathLength.push_back(in_pathLength + 5100*(randUniformReal(randGenerator)-0.5));
    vector<double> linDensity;
    cout << "Enter linDensities (0 to terminate): ";
    while ((cin >> input) && input != 0)
        linDensity.push_back(input);
    vector<double> persisLength;
    cout << "Enter persisLengths (0 to terminate): ";
    while ((cin >> input) && input != 0)
        persisLength.push_back(input);
    cout << "Enter segConvFactor: ";
    double segConvFactor;
    cin >> segConvFactor;
    WLCCollector myCollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             "data",
                             segConvFactor,
                             2.12,
                             false);
}