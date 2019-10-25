import numpy as np
import os
import logging
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis,\
    BrightestFluxWeightHypothesis
from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue
from alertstack.stats import GammaDistribution
from examples.fermi_neutrino_alert import blazar_analysis

if __name__ == "__main__":

    logging.getLogger().setLevel("INFO")

    blazar_analysis.iterate_run(
        n_trials=100,
        injection_hypo=AverageFluxWeightHypothesis,
        fraction=0.5,
        nsteps=5,
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


    for step, res in all_res.items():
        print("Fraction of neutrino alerts correlated to source: {0} \n".format(step))

        bkgs = dict()

        for key, val in res.items():
            print(key, np.mean(val), np.median(val), np.std(val))
            val = np.array(val)

            for name, thresh in levels:
                print(thresh[key])
                print("Fraction above {0}: {1}".format(
                    name, np.sum(val > thresh[key])/float(len(val))))