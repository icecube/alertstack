import numpy as np
import os
import logging
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis
from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue
from alertstack.stats import GammaDistribution
import unittest

blazar_cache = os.path.join(os.path.dirname(__file__), "test_cache/")

expected_levels = [
    ('Background Median', {'average_flux_weight': 6.302686494522481}),
    ('3 Sigma Discovery Potential', {'average_flux_weight': 17.877173064634526}),
    ('5 Sigma Discovery Potential', {'average_flux_weight': 29.756258163220224})
]

expected_res = [
    [0.5, 0.0016, 0.0],
    [0.908, 0.156, 0.002],
    [0.972, 0.478, 0.066],
    [0.998, 0.828, 0.24],
    [1.0, 0.922, 0.514],
    [1.0, 0.978, 0.718]
]

class TestFullChain(unittest.TestCase):

    def setUp(self):
        pass

    def test_chain(self):
        logging.info("Testing full analysis chain.")

        blazar_analysis = Analyse(
            Fermi4FGLBlazarCatalogue(),
            [AverageFluxWeightHypothesis],
            CircularisedNeutrinoAlertCatalogue(),
            cache_dir=blazar_cache,
            clean_cache=True
        )

        blazar_analysis.iterate_run(
            n_trials=100,
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
            bkg_median = np.median(val)
            bkg_expected = expected_levels[0][1][key]
            self.assertAlmostEqual(bkg_median/bkg_expected, 1.0, delta=0.05)
            sens_threshold[key] = bkg_expected

            gd = GammaDistribution(val)

            disc_3_res = gd.calculate_discovery_potential(3.)
            disc_3_expected =expected_levels[1][1][key]
            self.assertAlmostEqual(disc_3_res/disc_3_expected, 1.0, delta=0.1)
            disc_3_threshold[key] = disc_3_expected

            disc_5_res = gd.calculate_discovery_potential(5.)
            disc_5_expected = expected_levels[2][1][key]
            self.assertAlmostEqual(disc_5_res / disc_5_expected, 1.0, delta=0.1)
            disc_5_threshold[key] = disc_5_expected

        levels = [
            ("Background Median", sens_threshold),
            ("3 Sigma Discovery Potential", disc_3_threshold),
            ("5 Sigma Discovery Potential", disc_5_threshold)
        ]

        reduced_chi2 = 0.
        ndof = 0.

        for i, (step, res) in enumerate(all_res.items()):
            logging.debug(f"Fraction of ASTROPHYISCAL neutrino alerts correlated to source: {step} \n")

            for key, val in res.items():
                logging.debug(f"{key}, {np.mean(val)}, {np.median(val)}, {np.std(val)}")
                val = np.array(val)

                for j, (name, thresh) in enumerate(levels):
                    logging.debug(thresh[key])

                    frac = np.sum(val > thresh[key])/float(len(val))

                    n_obs = float(len(val))
                    n_det = np.sum(val > thresh[key])

                    exp = expected_res[i][j] * n_obs

                    chi2 = ((n_det - exp)**2.)/(exp + 1.)

                    print(chi2, exp, n_det, exp**2., n_obs)

                    logging.debug(f"Fraction above {name}: {frac}")

                    print(expected_res[i][j], frac, abs(expected_res[i][j] - frac) > 0.06)

                    reduced_chi2 += chi2
                    ndof += 1

        print(f"Final {reduced_chi2}, {ndof}, {reduced_chi2/ndof}")

        self.assertLess(reduced_chi2/ndof, 1.2)

if __name__ == '__main__':
    logging.getLogger().setLevel("DEBUG")
    unittest.main()