import PolymerCppCore
import numpy as np

def getCppWLC(pathLength=1000.0, persisLength=1.0):
    """Generate the trajectory of an infinitesimal wormlike chain.

    Parameters
    ----------
    pathLength : float
        The length of the chain in atomic units, e.g. atoms,
        molecules, or basepairs.
    persisLength : float
        The persistence length of the chain in units of distance.

    Returns
    -------
    rawChain : array_like
        A Mx3 ordered array representing the chain's trajectory, where
        M is the number of points. Each row represents the x-, y-, and
        z-coordinates of the points in the trajectory.

    Examples
    --------
    Generate a single WLC that contains 100 points and has a
    persistence length of 10 line segments.

    >>> chain = getCppWLC(pathLength=100, persisLength=10)
    >>> print(chain.shape)
    (100, 3)
    
    """
    rawChain = np.array(PolymerCppCore.getWLC(pathLength, persisLength))   
    rawChain = np.reshape(rawChain, (-1, 3))
    return rawChain

def getCppSAWLC(pathLength=1000.0, persisLength=1.0, linkDiameter=0.5):
    """Generate a trajectory of a self-avoiding wormlike chain.

    Parameters
    ----------
    pathLength : float
        The length of the chain in atomic units, e.g. atoms,
        molecules, or basepairs.
    persisLength : float
        The persistence length of the chain in units of distance.

    Returns
    -------
    rawChain : array_like
        A Mx3 ordered array representing the chain's trajectory, where
        M is the number of points. Each row represents the x-, y-, and
        z-coordinates of the points in the trajectory.

    Examples
    --------
    Generate a single self-avoiding wormlike chain that is 100 units
    long, has a persistence length of 10, and a diameter that is 0.5
    times the length of a single link.

    >>> chain = getCppSAWLC(pathLength=100, persisLength=10, linkDiameter=0.5)
    >>> print(chain.shape)
    (100, 3)

    """
    rawChain = np.array(PolymerCppCore.getSAWLC(pathLength, persisLength, linkDiameter))
    rawChain = np.reshape(rawChain, (-1,3))
    return rawChain

'''
# These 2 functions return the gyration radius, you can get multiple
# values at once, but not for different values of pathLength at the
# same time. If you need a range of values of pathLength, just call
# this function N number of times with numpaths = 1.
def getCppWLCradii(numPaths=1,
                   pathLength= 1000.0,
                   linDensity=1.0,
                   persisLength=1.0,
                   bumped= True,
                   locPrecision=0.0):
	radii = np.array(PolymerCppCore.getWLCrgs(numPaths,
                                                  pathLength,
                                                  linDensity,
                                                  persisLength,
                                                  int(bumped),
                                                  locPrecision))
	if not bumped:
		return radii
	else: #search for the delimiting value of -1
		i = 0
		while (radii[i]>=0.0):
			i+=1
		return radii[i+1:-1] #return all values after the negative one

def getCppSAWLCradii(numPaths=1,
                     pathLength=1000.0,
                     linDensity=1.0,
                     persisLength=1.0,
                     segConvFactor=1.0,
                     bumped=True,
                     locPrecision=0.0,
                     linkDiameter=0.5):
	radii = np.array(PolymerCppCore.getSAWLCrgs(numPaths,
                                                    pathLength,
                                                    linDensity,
                                                    persisLength,
                                                    segConvFactor,
                                                    int(bumped),
                                                    locPrecision,
                                                    linkDiameter))
	if not bumped:
		return radii
	else: #search for the delimiting value of -1
		i = 0
		while (radii[i]>=0.0):
			i+=1
		return radii[i+1:-1] #return 
'''
