import numpy as np
import os
import logging
import argparse
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis,\
    BrightestFluxWeightHypothesis
#from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue
from alertstack.stats import GammaDistribution
from examples.fermi_blazar_neutrino_alert import blazar_analysis

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate TS distributions')
    parser.add_argument('--n_trials', type=int, default=500,
                            help = 'Number of trials')
    parser.add_argument('--fraction', type=float, default=0.5,
                            help = 'Maximum fraction of neutrinos to be correlated')
    parser.add_argument('--n_steps', type=int, default=10, help ='Number of steps')
    args = parser.parse_args()

    logging.getLogger().setLevel("INFO")

    blazar_analysis.iterate_run(
        n_trials=args.n_trials,
        injection_hypo=AverageFluxWeightHypothesis,
        fraction=args.fraction,
        n_steps=args.n_steps,
    )

    all_res = blazar_analysis.load_results()

    sens_threshold = dict()
    disc_3_threshold = dict()
    disc_5_threshold = dict()

    zero_key = 0.0

    for key, val in all_res[zero_key].items():
        sens_threshold[key] = np.median(val)
        gd = GammaDistribution(val)
        disc_3_threshold[key] = gd.calculate_discovery_potential(3.)
        disc_5_threshold[key] = gd.calculate_discovery_potential(5.)
        # print(Chi2(val))
        # input("?")

    levels = [
        ("Background Median", sens_threshold),
        ("3 Sigma Discovery Potential", disc_3_threshold),
        ("5 Sigma Discovery Potential", disc_5_threshold)
    ]

    above = dict()

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

    f1 = interpolate.interp1d(fracs, sens, kind='cubic')
    f2 = interpolate.interp1d(fracs, sig3, kind='cubic')
    f3 = interpolate.interp1d(fracs, sig5, kind='cubic')

    tmp = [i.weight for i in blazar_analysis.fixed_sources]
    avg_signalness = np.mean(tmp)
    n_events = len(tmp)
    print("\n##### Calculate sensitivity and discovery potential with {0} neutrino alerts (average signalness: {1:.1f} %) ######\n".format(n_events, 100*avg_signalness))

    for x in np.arange(0, args.fraction, 0.001): 
        if abs(f1(x)-0.9)<0.01 :
            print("Sensitivity at {0} of flux, expectation of {1:.2f}/{2} = {3:.2f}".format(
                x,(x*avg_signalness)*n_events,n_events,(x*avg_signalness)))
        if abs(f2(x)-0.5)<0.01 :
            print("3 Sigma discovery at {0} of flux, expectation of {1:.2f}/{2} = {3:.2f}".format(
                x,(x*avg_signalness)*n_events,n_events,(x*avg_signalness)))
        if abs(f3(x)-0.5)<0.005 : 
            print("5 Sigma discovery at {0} of flux, expectation of {1:.2f}/{2} = {3:.2f}".format(
                x,(x*avg_signalness)*n_events,n_events,(x*avg_signalness)))