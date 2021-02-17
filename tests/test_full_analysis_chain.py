import numpy as np
import os
import logging
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis
from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue
from alertstack.stats import GammaDistribution

blazar_cache = os.path.join(os.path.dirname(__file__), "test_cache/")

blazar_analysis = Analyse(
    Fermi4FGLBlazarCatalogue(),
    [AverageFluxWeightHypothesis],
    CircularisedNeutrinoAlertCatalogue(),
    cache_dir=blazar_cache,
    clean_cache=True
)

blazar_analysis.iterate_run(
    n_trials=50,
    injection_hypo=AverageFluxWeightHypothesis,
    fraction=0.5,
    n_steps=5,
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

levels = [
    ("Background Median", sens_threshold),
    ("3 Sigma Discovery Potential", disc_3_threshold),
    ("5 Sigma Discovery Potential", disc_5_threshold)
]

for step, res in all_res.items():
    logging.debug(f"Fraction of ASTROPHYISCAL neutrino alerts correlated to source: {step} \n")

    bkgs = dict()

    for key, val in res.items():
        logging.debug(f"{key}, {np.mean(val)}, {np.median(val)}, {np.std(val)}")
        val = np.array(val)

        for name, thresh in levels:
            logging.debug(thresh[key])
            logging.debug(f"Fraction above {name}: {np.sum(val > thresh[key])/float(len(val))}")