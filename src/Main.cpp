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

    //* Set up simulation parameters
    Stopwatch swatch; swatch.set_mode(REAL_TIME);
    int numPaths = 150000;

    vector<double> pathLength(numPaths,11000);
    for (int i=0; i<pathLength.size(); i++)
    {
        pathLength[i] += 1500*(randNormalReal(randGenerator)-0.5);
        if (pathLength[i] < 7000)
            pathLength[i] = 7000;
        if (pathLength[i] > 15000)
            pathLength[i] = 15000;
    }
            
    vector<double> linDensity, persisLength;
    vector<double> linkDiameter {11.0};
    for (double i = 5.0; i<60.0; i+=5.0)
    {
        linDensity.push_back(i);
        persisLength.push_back(i);
    }
    double segConvFactor = 1.0/11.0;
    double locPrecision = 15.0;
    //*/


    // Run self-avoiding simulation
    swatch.start("sw");
    SACollector myCollector(numPaths,
                         pathLength,
                         linDensity,
                         persisLength,
                         linkDiameter,
                         "dataSA",
                         segConvFactor,
                         locPrecision,
                         false);
    myCollector.startCollector();

    // Run WLC simulation
    WLCCollector myWLCCollector(numPaths,
                         pathLength,
                         linDensity,
                         persisLength,
                         "dataWLC",
                         segConvFactor,
                         locPrecision,
                         false);
    myWLCCollector.startCollector();
    swatch.stop("sw");
    cout << "Total simulation time: "
         << swatch.get_total_time("sw") << std::endl;
    /* SWITCH: Test case 1: Test sphere sampling

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
    //*/

    /* SWITCH: Test case 2: test for vector normalization
    Eigen::Vector3d point(6.3,0.7,-5.6);
    cout << "Vector: " << point << endl;
    cout << "Squared norm: " << point.squaredNorm() << endl;
    point /= sqrt(point.squaredNorm());
    cout << "Normalized vector: " << point << endl;
    cout << "Squared norm: " << point.squaredNorm() << endl;
    //*/

    /* SWITCH: Test case 3: create a single random walk
    Eigen::Vector3d startPoint(0.0,1.0,0.0);
    WormlikeChain WLC(10.05, 25, &startPoint);
    WLC.makeNewPath();
    for (auto & point: WLC.path)
    {
        cout << point(0) << ", " << point(1) << ", " << point(2) << endl;
    }
    //*/

    /* SWITCH: Test case 4: Create a WLCCollector,
               see if computed Rg matches theory
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
    //*/

    /* SWITCH: Test case 5: Create a single SAWLC
    vector<double> pathLength(1000,10000.0);

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
    //*/

    /* SWITCH: Test case 6: Create a SACollector,
               see if computed Rg matches theory

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
    //*/

    /* SWITCH: Actual computation algorithm.
    int numPaths = 20000;
    vector<double> pathLength(numPaths, 12000.0);
    vector<double> linDensity {5,7.5,10,12.5,15,20,25,30,35,40,50};
        //{45.0, 50.0, 55.0};
    vector<double> persisLength {0,1,3,5,7.5,10,12.5,15,17.5,20,
                                 22.5,25,27.5,30};
                                //35.0, 40.0, 45.0, 50.0, 70.0, 90.0,
                                //{55.0, 60.0, 65.0, 75.0, 80.0, 85.0,
                                //    95.0, 100.0, 105.0, 110.0};
    vector<double> linkDiameter {0.8};
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
    wlcCollector.startCollector();
    //*/

    /* SWITCH: Comparing resulting gyration radii for different algorithms
    //           with dependence on persistence length
    int numPaths = 2000;
    vector<double> pathLength(numPaths, 150.0);
    vector<double> linDensity {1};
    vector<double> persisLength {0,0.01,0.03,0.05,0.1,0.15,0.2,0.3,0.4,0.5,
                                0.75,1,1.25,2.5,4,6,8,10};
    //vector<double> linkDiameter {0.95};
    double segConvFactor = 1.0;
    double locPrecision = 1.0;
    for (double lD=0.1; lD<1.0; lD+=0.4)
    {
        vector<double> linkDiameter {lD};
        cout << lD << std::endl;
        for (auto & i : persisLength)
            cout << i << " ";
        cout << std::endl;
        {
            SACollector_Rosenbluth myCollector(numPaths,
                                 pathLength,
                                 linDensity,
                                 persisLength,
                                 linkDiameter,
                                 "dataSA8",
                                 segConvFactor,
                                 locPrecision,
                                 false);
            myCollector.startCollector();
        }
        cout << std::endl;
        {
            SACollector mySACollector(numPaths,
                                 pathLength,
                                 linDensity,
                                 persisLength,
                                 linkDiameter,
                                 "dataSA8",
                                 segConvFactor,
                                 locPrecision,
                                 false);
            mySACollector.startCollector();
        }
        cout << std::endl;
            WLCCollector wlcCollector(numPaths,
                                 pathLength,
                                 linDensity,
                                 persisLength,
                                 "data8",
                                 segConvFactor,
                                 locPrecision,
                                 false);
            wlcCollector.startCollector();
        cout << std::endl;
        for (auto & chainptr: wlcCollector.myChains)
            cout << chainptr->getTheoreticalRg() << " ";
        cout << std::endl << std::endl;
    }
    //*/

    /* SWITCH: Profiling computation time in dependence on pathLength
    Stopwatch swatch; swatch.set_mode(REAL_TIME);
    int numPaths = 5000;
    vector<double> linDensity {1};
    vector<double> persisLength {10};
    vector<double> linkDiameter {0.5};   
    double segConvFactor = 1.0;
    double locPrecision = 1.0;

    vector<double> pathLengths = {3, 10, 100, 200, 500, 1000, 1500, 2000, 2500,
                                3000, 5000, 7000, 10000, 13000, 15000};
    for (auto & pL : pathLengths)
    {
        cout << pL << " ";
        swatch.start("sw");
        vector<double> pathLength(numPaths, pL);

        SACollector mySACollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             linkDiameter,
                             "dataSA8",
                             segConvFactor,
                             locPrecision,
                             false);
        mySACollector.startCollector();
        swatch.stop("sw");
        cout << swatch.get_total_time("sw") << std::endl;
    }
    //*/


    /* SWITCH: Profiling computation time in dependence on persisLength
    Stopwatch swatch; swatch.set_mode(REAL_TIME);
    int numPaths = 5000;
    vector<double> linDensity {1};
    vector<double> persisLengths = {0,0.01,0.03,0.05,0.1,0.15,0.2,0.3,0.4,0.5,
                                0.75,1,1.25,2.5,4,6};
    vector<double> linkDiameter {0.5};
    double segConvFactor = 1.0;
    double locPrecision = 1.0;
    vector<double> pathLength(numPaths,3000);
    for (auto & pL : persisLengths)
    {
        cout << pL << " ";
        swatch.start("sw");
        vector<double> persisLength = {pL};

        SACollector mySACollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             linkDiameter,
                             "dataSA8",
                             segConvFactor,
                             locPrecision,
                             false);
        mySACollector.startCollector();
        swatch.stop("sw");
        cout << swatch.get_total_time("sw") << " ";

        swatch.reset("sw");
        swatch.start("sw");
        WLCCollector wlcCollector(numPaths,
                             pathLength,
                             linDensity,
                             persisLength,
                             "data8",
                             segConvFactor,
                             locPrecision,
                             false);
        wlcCollector.startCollector();
        swatch.stop("sw");
        cout << swatch.get_total_time("sw") << std::endl;
        swatch.reset("sw");
    }
    //*/


}