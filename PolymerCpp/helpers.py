"""Python wrappers to the C++ wormlike chain generation code.

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
Switzerland, Laboratory of Experimental Biophysics, 2017
See the LICENSE.txt file for more details.

"""

import PolymerCppCore
import numpy as np

def getCppWLC(pathLength=1000.0, persisLength=1.0):
    """Generate the trajectory of an infinitesimal wormlike chain.

    Parameters
    ----------
    pathLength : float
        The length of the chain in atomic units, e.g. atoms,
        molecules, base pairs, or segments.
    persisLength : float
        The persistence length of the chain in units of chain
        segments.

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
    (101, 3)
    
    """
    rawChain = np.array(PolymerCppCore.getWLC(pathLength, persisLength))

    while np.isnan(rawChain).any():
        rawChain = np.array(PolymerCppCore.getWLC(pathLength, persisLength))

    rawChain = np.reshape(rawChain, (-1, 3))

    return rawChain

def getCppWLC2D(pathLength=1000, persisLength=1.0):
    """Create a infinitesimally-thin wormlike chain in two dimensions.

    Parameters
    ----------
    pathLength : int
        The length of the chain in atomic units, e.g. atoms,
        molecules, base pairs, or segments.
    persisLength : float
        The persistence length of the chain in units of chain
        segments.

    Returns
    -------
    rawChain : array_like

    Examples
    --------
    Generates a single 2D WLC with 50 segments (51 vertices) and a
    persistence length of 10 segments.

    >>> chain = getCppWLC2D(pathLength=50, persisLength=10)
    >>> print(chain.shape)
    (51, 2)
    
    """
    rawChain = np.array(PolymerCppCore.getWLC2D(
        pathLength, persisLength))

    while np.isnan(rawChain).any():
        rawChain = np.array(PolymerCppCore.getWLC2D(
            pathLength,
            persisLength))

    rawChain = np.reshape(rawChain, (-1, 2))

    return rawChain

def getCppSAWLC(pathLength=1000.0, persisLength=1.0, linkDiameter=0.5):
    """Generate a trajectory of a self-avoiding wormlike chain.

    Parameters
    ----------
    pathLength : float
        The length of the chain in atomic units, e.g. atoms,
        molecules base pairs, or segments.
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
    (101, 3)

    """
    # Prevents exit of the Python interpreter if pathLength < 3 in C++ code
    if pathLength < 3:
        pathLengthCpp = 3
    else:
        pathLengthCpp = pathLength

    rawChain = np.array(PolymerCppCore.getSAWLC(pathLengthCpp,
                                                persisLength,
                                                linkDiameter))

    while np.isnan(rawChain).any():
        # TODO: Fix NaN's in C++ code
        rawChain = np.array(PolymerCppCore.getSAWLC(pathLengthCpp,
                                                    persisLength,
                                                    linkDiameter))
    
    rawChain = np.reshape(rawChain, (-1,3))[0:pathLength + 1, :]

    return rawChain

def radius_of_gyration(chain):
    """Compute the radius of gyration of a chain.

    This function computes the radius of gyration of a set of
    points. The points are defined as a NumPy array where each row
    corresponds to the Cartesian coordinate values and each column is
    one of the principle directions.

    Parameters
    ----------
    chain : array_like
        A MxN array of numbers representing the Cartesian coordinates
        of the points in the set where M is the number of points and
        N is the number of dimensions.

    Returns
    -------
    Rg : float
        The radius of gyration of the chain.

    """
    secondMoments = np.var(chain, axis = 0)
    Rg = (np.sum(secondMoments)) ** (0.5)

    return Rg

def end_to_end_distance(chain):
    """Compute the end-to-end distance of a chain.

    Parameters
    ----------
    chain : array_like
        A MxN array of numbers representing the Cartesian coordinates
        of the points in the set where M is the number of points and
        N is the number of dimensions.

    Returns
    -------
    R : float
        The end-to-end distance of the chain.

    """
    R = np.linalg.norm(chain[-1,:] - chain[0,:])

    return R

def theory_Rg_WLC(contour_length, persistence_length, dim=3):
    """The mean squared radius of gyration of the wormlike chain.

    Note that this returns the square root of the mean squared radius
    of gyration.

    Parameters
    ----------
    contour_length : float
        The total contour length of the chain
    persistence_length: float
        The persistence length of the chain, which is related to the
        chain's stiffness.
    dim : int
        The dimension of the chain. This may be either 2 or 3.

    Returns
    -------
    Rg : float
        The square root of the mean squared radius of gyration.
    
    """
    if dim == 2:
        term1 = 2 * contour_length * persistence_length / 3
        term2 = 4 * persistence_length**2
        term3 = 16 * persistence_length**3 / contour_length**2
        term4 = 2 * persistence_length * (1 - np.exp(-contour_length / 2 / persistence_length))
    elif dim == 3:
        term1 = contour_length * persistence_length / 3
        term2 = persistence_length**2
        term3 = 2 * persistence_length**3 / contour_length**2
        term4 = persistence_length * (1 - np.exp(-contour_length / persistence_length))
        
    Rg = np.sqrt(term1 - term2 + term3 * (contour_length - term4))
    
    return Rg

def theory_R_WLC(contour_length, persistence_length, dim=3):
    """The mean squared end-to-end distance of the wormlike chain.

    Note that this returns the square root fo the mean squared end-to-
    end distance.

    Parameters
    ----------
    contour_length : float
        The total contour length of the chain.
    persistence_length: float
        The persistence length of the chain, which is related to the
        chain's stiffness.
    dim : int
        The dimension of the chain. This may be either 2 or 3.

    Returns
    -------
    R : float
        The square root of the mean squared end-to-end distance.

    """
    if dim == 2:
        term1 = 4 * persistence_length * contour_length
        term2 = 2 * persistence_length / contour_length
        term3 = 1 - np.exp(-contour_length / 2 / persistence_length)
    elif dim == 3:
        term1 = 2 * persistence_length * contour_length
        term2 = persistence_length / contour_length
        term3 = 1 - np.exp(-contour_length / persistence_length)
    else:
        raise ValueError('dim must be either 2 or 3. {} was supplied'.format(dim))
            
    R = np.sqrt(term1 * (1 - term2 * term3))

    return R
