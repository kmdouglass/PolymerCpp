"""Tools to ensure the chain generation algorithms work as expected.

© All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE,
Switzerland, Laboratory of Experimental Biophysics, 2017
See the LICENSE.txt file for more details.

"""

import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
from PolymerCpp.helpers import getCppWLC, radius_of_gyration
from PolymerCpp.helpers import theory_Rg_WLC, theory_R_WLC, end_to_end_distance

def verifyWLC(algorithm,
              numChains=10,
              numExperiments=1000,
              pathLength=1000,
              persisLength=50,
              **kwargs):
    """Compares outputs of a wormlike chain algorithm to theory.

    This function performs two numerical experiments that verify the
    accuracy of the simulations. In the first, it performs a number
    `numExperiments` of experiments with `numChains` realizations in
    each experiment. The mean end-to-end distance and mean radius of
    gyration is computed for each group and a histogram of the results
    is displayed. The x-axis is the difference between the computed
    mean and the theoretical mean.

    In the second experiment, a number of chains equal to `numChains`
    are simulated for a range of chain contour lengths.  The mean
    end-to-end distance and mean radius of gyration for each group is
    then plotted vs. the contour length, along with the theoretical
    predications. The error bars denote the standard deviation of the
    group.

    The mean value of all the end-to-end distances and radii of
    gyration are output to the console.

    Parameters
    ----------
    algorithm : func
        The algorithm for generating a wormlike chain.
    numChains : int
        The number of chains per numerical experiment to generate.
    numExperiments : int
        The number of experiments to perform for calculating the
        mean-squared moments.
    pathLength : float
        The length of each chain in units of segments.
    persisLength : float
        The persistence length of the wormlike chain in units of
        chain segments.
    kwargs : dict
        A list of additional keyword arguments to pass to the
        algorithm function.

    Notes
    -----
    If a self-avoiding wormlike chain is provided as the argument to
    `algorithm`, the results will be compared to the theory for the
    infinitesimally thin wormlike chain.

    """
    
    # Experiment 1: Compute mean squared radius of gyration
    mean_rgs1 = np.zeros(numExperiments)
    mean_rs1 = np.zeros(numExperiments)
    for i in range(numExperiments):

        rg = np.zeros(numChains)
        r  = np.zeros(numChains)

        # Generate the individual chains
        for j in range(numChains):

            rawChain = algorithm(pathLength,
                                 persisLength,
                                 **kwargs)
            rg[j] = radius_of_gyration(rawChain)
            r[j]  = end_to_end_distance(rawChain)

        mean_rgs1[i] = np.sqrt(np.mean(rg**2))
        mean_rs1[i] = np.sqrt(np.mean(r**2))

    # Experiment 2: Scaling behavior
    Lc = np.round(np.logspace(0, np.log10(pathLength)))
    mean_rgs2 = np.zeros(len(Lc))
    std_rgs2  = np.zeros(len(Lc))

    mean_rs2 = np.zeros(len(Lc))
    std_rs2  = np.zeros(len(Lc))
    for i in range(len(Lc)):

        rg = np.zeros(numChains)
        r  = np.zeros(numChains)
        
        for j in range(numChains):
            
            rawChain = algorithm(int(Lc[i]), persisLength, **kwargs)
            rg[j] = radius_of_gyration(rawChain)
            r[j]  = end_to_end_distance(rawChain)

        mean_rgs2[i] = np.sqrt(np.mean(rg**2))
        std_rgs2[i]  = np.std(rg)

        mean_rs2[i] = np.sqrt(np.mean(r**2))
        std_rs2[i]  = np.std(r)


    # Get the chain dimension
    dim = rawChain.shape[1]
        
    # Output results
    theory_rg = theory_Rg_WLC(pathLength, persisLength, dim=dim)
    theory_r = theory_R_WLC(pathLength, persisLength, dim=dim)

    print('Experimental inputs')
    print('-------------------')
    print('Number of chains:\t\t{0:d}\nNumber of experiments:\t\t{1:d}'.format(numChains, numExperiments))
    print('Contour length (for moments):\t{0:d}\npersistence length:\t\t{1:0.4f}'.format(pathLength, persisLength))
    print('Chain dimension:\t\t{0:d}'.format(dim))
    print('')
    print('Outputs')
    print('-------')
    print('Calculated <R>:\t\t\t{0:.4f}\nTheoretical WLC <R>:\t\t{1:.4f}'.format(np.mean(mean_rs1), theory_r))
    print('')
    print('Calculated <Rg>:\t\t{0:.4f}\nTheoretical WLC <Rg>:\t\t{1:.4f}'.format(np.mean(mean_rgs1), theory_rg))
    
    # Plot results
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(9, 5))
    
    ax0.hist(mean_rs1 - theory_r, bins=50, label=r'$ \langle R^2 \rangle^{1/2} $')
    ax0.hist(mean_rgs1 - theory_rg, bins=50, label=r'$ \langle R_g^2 \rangle^{1/2} $')
    ax0.set_title('Moments')
    ax0.set_xlabel(r'$ \Delta \langle R^2 \rangle^{1/2} , \, \Delta \langle R_g^2 \rangle^{1/2} } $')
    ax0.set_ylabel('Frequency')
    ax0.legend()

    ax1.errorbar(Lc, mean_rs2, yerr=std_rs2, fmt='o', label=r'$ \langle R^2 \rangle_{numeric} $', markersize=2)
    ax1.errorbar(Lc, mean_rgs2, yerr=std_rgs2, fmt='o', label=r'$ \langle R_g^2 \rangle_{numeric} $', markersize=2)
    ax1.plot(Lc, theory_R_WLC(Lc, persisLength, dim=dim), ':k', label=r'WLC $ \langle R^2 \rangle_{theory} $')
    ax1.plot(Lc, theory_Rg_WLC(Lc, persisLength, dim=dim), '--k', label=r'WLC $ \langle R_g^2 \rangle_{theory} $')
    ax1.set_title('Scaling behavior')
    ax1.set_xlabel('Contour length')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.legend()

    for ax in (ax0, ax1):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')
        
    plt.show()

if __name__ == '__main__':
    from PolymerCpp.helpers import getCppWLC
    verifyWLC(getCppWLC)
