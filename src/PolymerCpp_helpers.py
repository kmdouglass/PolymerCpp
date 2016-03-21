import PolymerCpp
import numpy as np

def getCppWLC(pathLength = 1000.0,\
			  linDensity = 1.0,\
			  persisLength = 1.0,\
			  bumped = False,\
			  locPrecision = 0.0):
	rawChain = np.array( PolymerCpp.getWLC(pathLength, linDensity,\
			    persisLength, int(bumped), locPrecision) )
	return np.reshape(rawChain, (-1,3))

def getCppSAWLC(pathLength = 1000.0,\
			  linDensity = 1.0,\
			  persisLength = 1.0,\
			  bumped = False,\
			  locPrecision = 0.0,\
			  linkDiameter = 0.5):
	rawChain = np.array( PolymerCpp.getSAWLC(pathLength, linDensity,\
			    persisLength, int(bumped), locPrecision,\
			    linkDiameter) )
	return np.reshape(rawChain, (-1,3))


def getCppWLCradii(numPaths = 10,\
				pathLength = 1000.0,\
			  linDensity = 1.0,\
			  persisLength = 1.0,\
			  bumped = True,\
			  locPrecision = 0.0):
	radii = np.array( PolymerCpp.getWLCrgs(numPaths, pathLength, linDensity,\
			    persisLength, int(bumped), locPrecision) )
	if not bumped:
		return radii
	else: #search for the delimiting value of -1
		i = 0
		while (radii[i]>=0.0):
			i+=1
		return radii[i+1:-1] #return all values after the negative one

def getCppSAWLCradii(numPaths = 10,\
				pathLength = 1000.0,\
			  linDensity = 1.0,\
			  persisLength = 1.0,\
			  bumped = True,\
			  locPrecision = 0.0,\
			  linkDiameter = 0.5):
	radii = np.array( PolymerCpp.getSAWLCrgs(numPaths, pathLength, linDensity,\
			    persisLength, int(bumped), locPrecision, linkDiameter) )
	if not bumped:
		return radii
	else: #search for the delimiting value of -1
		i = 0
		while (radii[i]>=0.0):
			i+=1
		return radii[i+1:-1] #return 