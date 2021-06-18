from astropy.io import fits
import numpy as np
import os
import logging
from alertstack import IsotropicExtragalacticCatalogue, Hypothesis, is_outside_GP, alertstack_data_dir
from numpy.lib.recfunctions import rename_fields

class Fermi4FGLBlazarCatalogue(IsotropicExtragalacticCatalogue):

    @staticmethod
    def parse_data():

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        with fits.open(os.path.join(alertstack_data_dir, "table_4LAC.fits")) as hdul:
            cat = hdul[1].data
        cat = np.sort(cat, order="Flux1000")[::-1]

        logging.info("Selecting blazars from 4FGL catalogue")

        blazar_class = ["bll", "BLL", "fsrq", "FSRQ", "bcu", "BCU"]

        logging.info("Using all sources from class {0}".format(blazar_class))
        #
        mask = np.array([x["CLASS"] in blazar_class for x in cat])
        blazars = np.array(cat[mask])

        cut_e = -11.6
        mask_e = np.array([x["Energy_Flux100"]>10**cut_e for x in blazars])
        blazars = np.array(blazars[mask_e])
        
        maps = [
            ("RAJ2000", "ra_rad"),
            ("DEJ2000", "dec_rad"),
        ]

        for (old_key, new_key) in maps:

            blazars = rename_fields(blazars, {old_key: new_key})

        mask_GP = [is_outside_GP(blazars['ra_rad'][i],blazars['dec_rad'][i]) for i in range(len(blazars))]
        blazars = blazars[mask_GP] 

        logging.info("Found {0} sources in total".format(len(blazars)))

        return blazars

    def scramble(self):
        ra, dec = self.scramble_positions()
        cat = np.copy(self.data)
        cat['ra_rad'] = ra
        cat["dec_rad"] = dec
        return cat


class AverageFluxWeightHypothesis(Hypothesis):
    name = "average_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        return cat_data["Flux1000"]

class BrightestFluxWeightHypothesis(Hypothesis):
    name = "brightest_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        weights = cat_data["Flux1000"]
        weights[100:] = 0.
        return weights

