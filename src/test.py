import PolymerCpp as PCPP
from PolymerCpp_helpers import *
import numpy as np
import time
import matplotlib.pyplot as plt

lens = (1235.0, 354.0, 3549.0, 7567.3)

a = getCppSAWLC()
then = time.time()
a = getCppSAWLCradii(numPaths = 10000, pathLength = 1000, persisLength = 0.1, linkDiameter = 0.75)
now = time.time()
print("Time elapsed 1: "+str(now-then))

plt.hist(a, bins = np.arange(0,100,0.5))

a = np.array([])
then = time.time()
for i in range(10000):
	a = getCppSAWLCradii(numPaths = 1, pathLength = 1000, persisLength = 0.1, linkDiameter = 0.75)
now = time.time()
print("Time elapsed 2: "+str(now-then))

plt.show()
