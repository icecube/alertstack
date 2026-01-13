import numpy as np
import os
import logging
import argparse
from scipy import interpolate
from scipy.optimize import bisect
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis
from alertstack.stats import GammaDistribution
from examples.fermi_blazar_neutrino_alert import blazar_analysis, plot_ts

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate TS distributions')
    parser.add_argument('--n_trials', type=int, default=500,
                            help = 'Number of trials')
    parser.add_argument('--fraction', type=float, default=0.5,
                            help = 'Maximum fraction of neutrinos to be correlated')
    parser.add_argument('--n_steps', type=int, default=10, help ='Number of steps')
    args = parser.parse_args()
    
    '''
    n_trials: Number of trials to run for each injection strength. 10x this number
    will be run as background trials
    fraction: Maximum fraction of astrophysical neutrinos to be injected
    n_steps: Number of different injection steps to test, between 0 and fraction.
    '''

    logging.getLogger().setLevel("INFO")

    # Run analysis and save results
    blazar_analysis.iterate_run(
        n_trials=args.n_trials,
        injection_hypo=AverageFluxWeightHypothesis,
        fraction=args.fraction,
        n_steps=args.n_steps,
    )

    # Load the latest results (search in ./cache) 
    all_res = blazar_analysis.load_results()

    sens_threshold = dict()
    disc_3_threshold = dict()
    disc_5_threshold = dict()

    zero_key = 0.0

    # Calculate threshold for sensitivity and discovery potential by fitting the 
    # background TS distribution (no correlation injected) to a gamma distribution
    for key, val in all_res[zero_key].items():
        sens_threshold[key] = np.median(val)
        gd = GammaDistribution(val)
        disc_3_threshold[key] = gd.calculate_discovery_potential(3.)
        disc_5_threshold[key] = gd.calculate_discovery_potential(5.)
        
        # Plot TS Distribution
        plot_ts(val, sens_threshold[key], disc_3_threshold[key], disc_5_threshold[key], gd)

    levels = [
        ("Background Median", sens_threshold),
        ("3 Sigma Discovery Potential", disc_3_threshold),
        ("5 Sigma Discovery Potential", disc_5_threshold)
    ]

    above = dict()
    
    # Calculate for each step the fraction of trials that are above a certain threshold
    # and store the information in the 'above' dictionary
    for step, res in all_res.items():
        print("\nFraction of neutrino alerts correlated to source: {0} \n".format(step))

        bkgs = dict()
        temp = []

        for key, val in res.items():
            print(key, np.mean(val), np.median(val), np.std(val))
            val = np.array(val)

            for name, thresh in levels:
                print(thresh[key])
                print("Fraction above {0}: {1}".format(
                    name, np.sum(val > thresh[key])/float(len(val))))
                temp.append(np.sum(val > thresh[key])/float(len(val)))        
        above[step] = temp

    fracs = list(above.keys())
    sens = [list(above.values())[i][0] for i in range(len(fracs))]
    sig3 = [list(above.values())[i][1] for i in range(len(fracs))]
    sig5 = [list(above.values())[i][2] for i in range(len(fracs))]
    
    # Interpolate a curve to the data points of the fraction of trials above each threshold
    f1 = interpolate.interp1d(fracs, sens, kind='cubic')
    f2 = interpolate.interp1d(fracs, sig3, kind='cubic')
    f3 = interpolate.interp1d(fracs, sig5, kind='cubic')
    print(f"\n\n\n*\t*\t*\t*\t*\t*\t*\n\nsens: {sens}\nsig3: {sig3}\nsig5: {sig5}\n\n*\t*\t*\t*\t*\t*\t*\n\n\n")

    # Calculate average signalness and number of neutrino alerts
    tmp = [i.weight for i in blazar_analysis.fixed_sources]
    avg_signalness = np.mean(tmp)
    n_events = len(tmp)

    # Calculate flux needed to achieve sensitivity and discovery potential 
    print("\n------- Sensitivity and discovery potential with {0} neutrino alerts (average signalness: {1:.1f} %) --------\n".format(n_events, 100*avg_signalness))

    x1 = bisect(lambda x: f1(x)-0.9, 0, args.fraction, xtol=1e-6)
    x2 = bisect(lambda x: f2(x)-0.5, 0, args.fraction, xtol=1e-6)
    x3 = bisect(lambda x: f3(x)-0.5, 0, args.fraction, xtol=1e-6)

    print("Sensitivity at {0:.3f} of flux, expectation of {1:.1f}/{2} = {3:.2f}".format(
                x1,(x1*avg_signalness)*n_events,n_events,(x1*avg_signalness)))
    print("3 Sigma discovery at {0:.3f} of flux, expectation of {1:.1f}/{2} = {3:.2f}".format(
                x2,(x2*avg_signalness)*n_events,n_events,(x2*avg_signalness)))
    print("5 Sigma discovery at {0:.3f} of flux, expectation of {1:.1f}/{2} = {3:.2f}".format(
                x3,(x3*avg_signalness)*n_events,n_events,(x3*avg_signalness)))
    print("\n----------------------------------------------------------------------------------------------------------\n")
