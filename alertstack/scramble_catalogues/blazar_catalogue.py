from astropy.io import fits
import numpy as np
import os
import logging
from pathlib import Path
import pickle as pkl
from alertstack import IsotropicExtragalacticCatalogue, Hypothesis, is_outside_GP, alertstack_data_dir
from numpy.lib.recfunctions import rename_fields
from astropy import units as u

# For the LC 

try:
    a = Path(alertstack_data_dir) / 'weights_LC.pkl'
    with a.open('rb') as f:
        data_lc = pkl.load(f)
except KeyError:
    logging.warning("The weights from the light curves could not be loaded. If you do not have them, importing MonthlyFluxWeightHypothesis will raise an error.")


class Fermi4FGLBlazarCatalogue(IsotropicExtragalacticCatalogue):
    '''
    Loads Fermi 4LAC-DR, selects blazars and applies a cut on the energy flux and on the latitude.
    '''

    @staticmethod
    def parse_data():

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        # Load catalog
        with fits.open(os.path.join(alertstack_data_dir, "table-4LAC-DR2-h.fits")) as hdul:
            cat = hdul[1].data
        cat = np.sort(cat, order="Energy_Flux100")[::-1]

        # Select blazars
        logging.info("Selecting blazars from 4FGL catalogue")

        blazar_class = ["bll", "BLL", "fsrq", "FSRQ", "bcu", "BCU"]

        logging.info("Using all sources from class {0}".format(blazar_class))
        mask = np.array([x["CLASS"] in blazar_class for x in cat])
        blazars = np.array(cat[mask])

        # Apply cut on energy flux
        cut_e = -11.6
        mask_e = np.array([x["Energy_Flux100"]>10**cut_e for x in blazars])
        blazars = np.array(blazars[mask_e])

        maps = [
            ("RAJ2000", "ra_rad"),
            ("DEJ2000", "dec_rad"),
        ]

        for (old_key, new_key) in maps:

            blazars = rename_fields(blazars, {old_key: new_key})

        # Apply cut on latitude 
        mask_GP = [is_outside_GP(blazars['ra_rad'][i],blazars['dec_rad'][i]) for i in range(len(blazars))]
        blazars = blazars[mask_GP] 

        logging.info("Found {0} sources in total".format(len(blazars)))

        return blazars 

    def scramble(self):
        ra, dec = self.scramble_positions_outside_GP()
        cat = np.copy(self.data)
        cat['ra_rad'] = ra
        cat["dec_rad"] = dec
        return cat
    


class AverageFluxWeightHypothesis(Hypothesis):
    name = "average_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        try:
            return cat_data["Energy_Flux100"]
        except:
            return cat_data['X band map']

        
class BrightestFluxWeightHypothesis(Hypothesis):
    name = "brightest_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        weights = cat_data["Energy_Flux100"]
        weights[100:] = 0.
        return weights
    
    
class MonthlyFluxWeightHypothesis(Hypothesis):
    '''
    In this hypothesis the blazar is weighted by the flux in the monthly time 
    bin of the neutrino arrival time (nu_at).
    '''
    name = "monthly_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data, nu_at, ignore_times=False):

        weights = [data_lc[name][nu_at] for name in cat_data['Source_Name']]
                
        return weights
    