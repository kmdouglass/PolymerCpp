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
    /*
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
    
*/
    // Test case 5: Create a single SAWLC
    /*vector<double> pathLength(1000,10000.0);

    Eigen::Vector3d startPoint(0.0,1.0,0.0);
    SAWLC selfAvoidingWLC(1000,
                          pathLength,
                          1,
                          20,
                          1.0/10,
                          2.1,
                          2.5,
                          &startPoint);
    selfAvoidingWLC.makePath(pathLength[0]);
    for (auto & point: selfAvoidingWLC.path)
    {
        cout << point(0) << ", " << point(1) << ", " << point(2) << endl;
    }
    */

    // Test case 6: Create a SACollector,
    //              see if computed Rg matches theory
    /*
    int numPaths; double input;
    cout << "Enter numPaths: ";
    cin >> numPaths;
    vector<double> pathLength;
    double in_pathLength;
    cout << "Enter base pathLength: ";
    cin >> in_pathLength;
    double in_pathLengthSpread;
    cout << "Enter pathLength spread: ";
    cin >> in_pathLengthSpread;
    for (int i=0; i<numPaths; i++)
        pathLength.push_back(in_pathLength + in_pathLengthSpread *
                    (randUniformReal(randGenerator)-0.5));
    vector<double> linDensity;
    cout << "Enter linDensities (0 to terminate): ";
    while ((cin >> input) && input != 0)
        linDensity.push_back(input);
    vector<double> persisLength;
    cout << "Enter persisLengths (0 to terminate): ";
    while ((cin >> input) && input != 0)
        persisLength.push_back(input);
    vector<double> linkDiameter;
    cout << "Enter linkDiameters (0 to terminate): ";
    while ((cin >> input) && input != 0)
        linkDiameter.push_back(input);
    cout << "Enter segConvFactor: ";
    double segConvFactor;
    cin >> segConvFactor;
    SACollector myCollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             linkDiameter,
                             "dataSA",
                             segConvFactor,
                             2.12,
                             false);
    myCollector.startCollector();
*/
    // Generating simulated data
    /*int numPaths = 20000;
    vector<double> pathLength(20000, 12000.0);
    vector<double> linDensity {5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.5, 
                            15.0, 17.5, 20.0, 22.5, 25.0, 30.0, 35.0, 40.0};
        //{45.0, 50.0, 55.0};
    vector<double> persisLength {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                                 10.0, 11.0, 13.0, 15.0, 17.5, 20.0, 25.0, 30.0};
                                //35.0, 40.0, 45.0, 50.0, 70.0, 90.0,
                                //{55.0, 60.0, 65.0, 75.0, 80.0, 85.0,
                                //    95.0, 100.0, 105.0, 110.0};
    vector<double> linkDiameter {0.02};
    double segConvFactor = 1.0;
    double locPrecision = 2.45;
    */
    int numPaths = 5;
    vector<double> pathLength(numPaths, 30000.0);
    vector<double> linDensity {15.0, 30, 49};
        //{45.0, 50.0, 55.0};
    vector<double> persisLength {0, 10, 100, 150};
                                //35.0, 40.0, 45.0, 50.0, 70.0, 90.0,
                                //{55.0, 60.0, 65.0, 75.0, 80.0, 85.0,
                                //    95.0, 100.0, 105.0, 110.0};
    vector<double> linkDiameter {0.02};
    double segConvFactor = 1.0;
    double locPrecision = 2.45;



    SACollector myCollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             linkDiameter,
                             "dataSA5",
                             segConvFactor,
                             locPrecision,
                             false);
    myCollector.startCollector();


    WLCCollector wlcCollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             "data5",
                             segConvFactor,
                             locPrecision,
                             false);
}