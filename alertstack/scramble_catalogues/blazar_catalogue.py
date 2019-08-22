from astropy.io import fits
import numpy as np
import logging
from alertstack import IsotropicExtragalacticCatalogue, Hypothesis
from numpy.lib.recfunctions import rename_fields

class Fermi4FGLBlazarCatalogue(IsotropicExtragalacticCatalogue):

    @staticmethod
    def parse_data():

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        hdul = fits.open("data/gll_psc_v19.fit")
        cat = hdul["LAT_Point_Source_Catalog"].data
        cat = np.sort(cat, order="Flux1000")[::-1]

        logging.info("Selecting blazars from 4FGL catalogue")

        blazar_class = ["bll", "BLL", "fsrq", "FSRQ"]

        logging.info("Using all sources from class {0}".format(blazar_class))

        mask = np.array([x["CLASS1"] in blazar_class for x in cat])
        blazars = cat[mask]

        maps = [
            ("RAJ2000", "ra_rad"),
            ("DEJ2000", "dec_rad"),
        ]

        for (old_key, new_key) in maps:

            blazars = rename_fields(blazars, {old_key: new_key})

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

