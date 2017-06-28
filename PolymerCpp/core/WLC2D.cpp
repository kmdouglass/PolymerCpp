#include "WLC2D.h"

using namespace std;

extern std::minstd_rand randGenerator;
extern std::normal_distribution<double> randNormalReal;


WLC2D::WLC2D(int in_pathLength,
             double in_persisLength, 
             Eigen::Vector2d * in_initPoint) 
    : Path2D(in_pathLength, in_initPoint)
{
    persisLength = in_persisLength;
}

void WLC2D::makePath(int in_pathLength)
{
    
    // Add one because two vertices = one segment
    int numVerts = in_pathLength + 1;
    // check if numVerts is in valid range
    if (numVerts < 1)
    {
        std::stringstream buffer;
        buffer << "The number of segments must be greater than 3, but a "
            << "value of " << numVerts << "was supplied. "
            << "Why would you want such a short chain anyways?" << std::endl;
        throw std::out_of_range(buffer.str());
    }

    Eigen::Vector2d currPoint = initPoint;
    Eigen::Vector2d nextPoint;

    // Define the set of rotation matrices that define the walk.
    double angle;
    vector<Eigen::Rotation2D<double>> rotations;

    // Since we don't have to check for collisions,
    // we can pre-generate all required rotations
    for (int i=0; i<numVerts; i++)
    {
        angle = pow(1 / persisLength, 0.5) * randNormalReal(randGenerator);
        Eigen::Rotation2D<double> rot2(angle);
        rotations.push_back(rot2);
    }

    // Primary iterative loop for creating the chain
    Eigen::Vector2d * workingPath = new Eigen::Vector2d[numVerts];
    workingPath[0] = Eigen::Vector2d(0.0, 0.0);
    workingPath[1] = initPoint;
    for (int i=2; i<numVerts; i++)
    {
        // Rotates the current position to the next
        nextPoint = rotations[i] * currPoint;

        // Append nextPoint to array of vectors on the path
        workingPath[i] = nextPoint;
        currPoint = nextPoint;
    }

    // Add up the vectors in path to create the polymer
    path[0] = workingPath[0];
    for (int i=1; i<numVerts; i++)
    {
        path.push_back(workingPath[i] + path[i-1]);
    }
    
    delete workingPath;
}