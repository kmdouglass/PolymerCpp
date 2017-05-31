import PolymerCppCore
import numpy as np

# These 2 functions return a chain with specified parameters in simulation units.
def getCppWLC(pathLength=1000.0,
              linDensity=1.0,
              persisLength=1.0,
              segConvFactor=1.0,
              bumped=False,
	      locPrecision=0.0):
	rawChain = np.array(PolymerCppCore.getWLC(pathLength,
                                                  linDensity,
                                                  persisLength,
                                                  segConvFactor,
                                                  int(bumped),
                                                  locPrecision))
	return np.reshape(rawChain, (-1,3))

def getCppSAWLC(pathLength=1000.0,
                linDensity=1.0,
                persisLength=1.0,
                segConvFactor=1.0,
		bumped=False,
		locPrecision=0.0,
		linkDiameter=0.5):
	rawChain = np.array(PolymerCppCore.getSAWLC(pathLength,
                                                    linDensity,
                                                    persisLength,
                                                    segConvFactor,
                                                    int(bumped),
                                                    locPrecision,
                                                    linkDiameter))
	return np.reshape(rawChain, (-1,3))


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
