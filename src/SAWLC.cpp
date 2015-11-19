#include "SAWLC.h"

using namespace std;

extern std::default_random_engine randGenerator;
extern std::uniform_real_distribution<double> randUniformReal;
extern std::normal_distribution<double> randNormalReal;

SAWLC::SAWLC(int in_numPaths, vector<double> & in_pathLength, 
                  double in_linDensity, double in_persisLength,
                  double in_linkDiameter, double in_segConvFactor, 
                  double in_locPrecision, Eigen::Vector3d * in_initPoint)
    : Path(in_numPaths, in_pathLength, in_linDensity,
            in_segConvFactor, in_initPoint)
{
	if (1.0/in_linDensity < in_linkDiameter)
    {
        std::stringstream buffer;
        buffer << "The link diameter is larger than the distance between links. "
            << std::endl << "Link diameter: " << in_linkDiameter << std::endl
            << "Distance between links: " << 1.0/in_linDensity << std::endl;
        throw std::out_of_range(buffer.str());
    }
    persisLength = in_persisLength;
    locPrecision = in_locPrecision;
    linkDiameter = in_linkDiameter;
    simLinkDiameter = in_linkDiameter * in_linDensity;
    simPersisLength = in_persisLength;
    defaultWeight = 1.0;
}

void SAWLC::makePath(double in_pathLength)
{   
    /* this chunk of code is the same in SAWLC and WLC
     * but I cant put it to Path:: because it would go out of scope
     */
    int numSegments = (int) (in_pathLength / linDensity);
    // check if numSegments is in valid range
    if (numSegments < 1)
    {
        std::stringstream buffer;
        buffer << "The number of segments must be greater than 1, but a "
            << "value of " << numSegments << "was supplied." << std::endl;
        throw std::out_of_range(buffer.str());
    }

    // split numSegments into integer and leftover
    // no leftover for the time being
    int numSeg;
    numSeg = (int)numSegments;

    double sigma = abs(simPersisLength)>0.00001 ? pow(2.0 / simPersisLength, 0.5) 
                                             : 999999999.9;
    Eigen::Vector3d currPoint = initPoint;
    Eigen::Vector3d dispVector, nextPoint;
    double projDistance;

    // Primary iterative loop for creating the chain
    vector<Eigen::Vector3d> workingPath(numSeg, Eigen::Vector3d::Constant(0.0));
    vector<Eigen::Vector3d> cumulativePath(numSeg, Eigen::Vector3d::Constant(0.0));
    workingPath[1] = initPoint;
    cumulativePath[1] = initPoint;
    vector<double> stepWeights(numSeg,0.0);
    stepWeights[0] = 1.0; stepWeights[1] = defaultWeight;
    vector<int> collisionPositions;
    double angDisp, tanPlaneDisp;
    Eigen::Vector3d * randVec = NULL;
    for (int i=2; i<numSeg; i++)
    {
    	// compute weight of step by integration over space
        getPossibleCollisions(cumulativePath, collisionPositions, i);

        stepWeights[i] = getNextStepWeight(cumulativePath, collisionPositions, i);
        // the second-to-last link is always available for collision,
        // we ignored it above, we add it now so we dont do it in the loop 
        collisionPositions.push_back(i-2);
        int j = 0;
        while (true)
        // keep generating vectors until you get one which is ok
        {
            j++;
            angDisp = sigma * randNormalReal(randGenerator);
            tanPlaneDisp = sin(angDisp);
            randVec = randPointSphere();
            dispVector = randVec->cross(currPoint);
            dispVector /= sqrt(dispVector.squaredNorm());
            dispVector = dispVector * tanPlaneDisp;
            projDistance = 1 - cos(angDisp);
            nextPoint = ((1 - projDistance) * currPoint) + dispVector;
            if (!checkCollision(nextPoint, cumulativePath, collisionPositions, i))
            {
                workingPath[i] = nextPoint;;
                cumulativePath[i] = nextPoint + cumulativePath[i-1];
                currPoint = nextPoint;
                delete randVec;
                break;
            }
            delete randVec;
            /*if (j%1000 == 0)
            {
                cout << "Tried new vector 1000 times, something is wrong!"
                    << std::endl; cout.flush(); 
            }*/
            if (j%10000 == 0)
            {
                cout << "DEAD END" << std::endl; cout.flush();
                throw 42;
            }
        }
    	// move to next step
    }

    // copy final vector to path
    path.clear();
    for (int i=0; i<numSeg; i++)
    {
        path.push_back(cumulativePath[i]);
    }

}

double SAWLC::getNextStepWeight(vector<Eigen::Vector3d> & cumulativePath,
                                vector<int> & collisionPositions,
                                int where)
{
    // check for possible collisions from other links

    if (collisionPositions.size() == 0)
        return defaultWeight;
    else
        return computeIntegral(cumulativePath, collisionPositions, where);
}

void SAWLC::getPossibleCollisions(
                    vector<Eigen::Vector3d> & cumulativePath,
                    vector<int> & collisionPositions,
                    int where)
{
    collisionPositions.clear();
    Eigen::Vector3d distanceVector; double distance;
    Eigen::Vector3d endPoint = cumulativePath.at(where-1);

    
    for (int i = 0; i < where-2; i++)
    // ignore last 2 points because they are always colliding
    {
        distanceVector = endPoint - cumulativePath[i];
        distance = sqrt(distanceVector.squaredNorm());
        if (distance < (1.0 + simLinkDiameter))
        {
            collisionPositions.push_back(i);
            /*cout << "Possible collision found at link no. "
                << where << ", with link no. "
                << i << std::endl; cout.flush();*/
        }
    }
    
    
    /*if (collisionPositions.size() > 0)
    {
        cout << "Found " << collisionPositions.size()
            << " possible collision(s) with body of chain! "
            << "At point " << cumulativePath.size()-1 << endl
            << collisionPositions.at(0) << endl; cout.flush();
    }*/
}

bool SAWLC::checkCollision(Eigen::Vector3d nextPoint,
                           vector<Eigen::Vector3d> & cumulativePath,
                           vector<int> & collisionPositions,
                           int where)
{
    Eigen::Vector3d endPoint = nextPoint + cumulativePath[where-1];
    //cout << "endPoint: " << std::endl << endPoint << std::endl;
    Eigen::Vector3d distanceVector; double distance;
    //cout << "collisionPositions size: " << collisionPositions.size() << std::endl;
    for (auto & position : collisionPositions)
    {
        distanceVector = endPoint - cumulativePath[position];
        distance = sqrt(distanceVector.squaredNorm());
        if (distance < (simLinkDiameter))
        {
            /*cout << "Collision averted at link no. "
                << where-1 << std::endl
                << "Distance: " << distance << std::endl
                << "simLinkDiameter: " << simLinkDiameter << std::endl; cout.flush();
            */
            return true;
        }
    }
    return false;
}

double SAWLC::computeIntegral(vector<Eigen::Vector3d> & cumulativePath,
                       vector<int> & collisionPositions, int where)
{
    return 1.0;
}



SACollector::SACollector(int in_numPaths,
                 vector<double> & in_pathLength,
                 vector<double> & in_linDensity,
                 vector<double> & in_persisLength,
                 vector<double> & in_linkDiameter,
                 string in_nameDB,
                 double in_segConvFactor,
                 double in_locPrecision,
                 bool in_fullSpecParam)
    : Collector(in_numPaths, in_pathLength, in_linDensity, in_persisLength,
                in_nameDB, in_segConvFactor, in_locPrecision, in_fullSpecParam)
{
    linkDiameter = in_linkDiameter;
    convSegments(linkDiameter, in_linkDiameter, segConvFactor, true);
    //copies variables to class
    if (!fullSpecParam)
    {
        vector<double> permLinkDiameter;
        vector<double> permLinDensity;
        vector<double> permPersisLength;
        for (int i=0; i<linDensity.size(); i++) {
            for (int j=0; j<linkDiameter.size(); j++) {
                permLinDensity.push_back(linDensity[i]);
                permPersisLength.push_back(persisLength[i]);
                permLinkDiameter.push_back(linkDiameter[j]);
            }
        }
        linDensity = permLinDensity;
        persisLength = permPersisLength;
        linkDiameter = permLinkDiameter;
    }
    for (int i=0; i<linDensity.size(); i++)
    {
        Eigen::Vector3d startDir(1.0,0.0,0.0);
        SAWLC * myChain = 
            new SAWLC(numPaths, pathLength, linDensity[i], 
                    persisLength[i], linkDiameter[i], segConvFactor,
                    locPrecision, &startDir);

        myChains.push_back(myChain);
    }

}

SAWLC * SACollector::getChainPointer(int i)
{
    return myChains[i];
}

SACollector::~SACollector()
{
    for (auto & chain: myChains)
    {
        delete chain;
    }
}

/*void SACollector::startCollector()
{
    cout << "Starting collector!" << endl;
    // Begin collecting wormlike chain conformation statistics.

    vector<RgDict*> data;
    #pragma omp parallel for schedule(static,1) num_threads(4)
         // parallelize the next loop
        for (int i=0; i<myChains.size(); i++) 
        {
            cout << "Simulation " << i << endl; cout.flush();
            data.push_back(parSimChain(myChains[i]));
        }
        // SAVE THE CALCULATED RgData
    // Opens or creates the database file for adding simulation results

    // sets maximum output precision
    fileDB << std::setprecision(std::numeric_limits<double>::digits10 + 1);
    fileDB.open(nameDB + ".txt", ios::out | ios::trunc);
    if (!fileDB.is_open())
    {
        std::stringstream buffer;
        buffer << "Could not open file: " << nameDB << ".txt" << std::endl;
        throw std::runtime_error(buffer.str());
    }


    // adds info about number of RgData objects to be saved
    fileDB << data.size();
    fileDB << " " << segConvFactor; // info about conversion factor

    for (auto & dict: data)
    {
        dict->addToDBfileFull(fileDB); // adds each object
        cout << "Density: " << dict->linDensity << ", Persistence length: "
             << dict->persisLength << std::endl;

        double sumRg = 0, sumRgBump = 0;
        for (int i=0; i<dict->Rg.size(); i++) {
            sumRg += dict->Rg.at(i);
            sumRgBump += dict->RgBump.at(i);
        }
        sumRg /= dict->Rg.size(); sumRgBump /= dict->RgBump.size();
        cout << "   Mean Rg:     " << sumRg << endl;
        cout << "   Mean RgBump: " << sumRgBump << endl;
        cout << "   Calculated:  " << theoreticalWLCRg(dict->linDensity,
                         dict->persisLength, dict->pathLength) << endl;
    }
    fileDB.close();


    // collect garbage
    for (auto & chain: myChains) delete chain;
}*/