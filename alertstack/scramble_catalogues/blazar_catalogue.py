import copy
import healpy as hp
import logging
import numpy as np
import pandas as pd
import pickle as pkl
import os

from alertstack import (
    alertstack_data_dir,
    Hypothesis,
    is_outside_GP,
    IsotropicExtragalacticCatalogue,
)
from astropy import units as u
from astropy.io import fits
from numpy.lib.recfunctions import rename_fields
from pathlib import Path

# Load the lightcurves (LC) to evaluate the monthly gamma-ray flux.
try:
    a = Path(alertstack_data_dir) / 'weights_LC.pkl'
    with a.open('rb') as f:
        data_lc = pkl.load(f)
except KeyError:
    logging.warning("The weights from the light curves could not be loaded. If you do not have them, importing MonthlyFluxWeightHypothesis will raise an error.")


class Fermi4FGLBlazarCatalogue(IsotropicExtragalacticCatalogue):
    ''' Loads Fermi 4LAC-DR, selects blazars, and applies a cut
    on the energy flux and on the latitude.
    '''

    @staticmethod
    def parse_data():
        ''' Loads Fermi 4LAC-DR, selects blazars, and applies a cut
        on the energy flux and on the latitude.
        '''

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        # Load catalog
        with fits.open(os.path.join(alertstack_data_dir, "table-4LAC-DR3-h.fits")) as hdul:
            cat = pd.DataFrame(hdul[1].data)
        for key in cat.keys():
            if cat[key].dtype == '>f8':
                cat[key] = cat[key].astype('f8')
        cat["Energy_Flux100"] = cat["Energy_Flux100"].astype('f8')
        cat = cat.sort_values("Energy_Flux100", ascending=False)

        # Select blazars
        logging.info("Selecting blazars from 4FGL catalogue")

        blazar_class = ["bll", "BLL", "fsrq", "FSRQ", "bcu", "BCU"]

        logging.info("Using all sources from class {0}".format(blazar_class))
        mask = np.array([df_class in blazar_class for df_class in cat["CLASS"]])
        blazars = cat[mask]

        # Apply cut on energy flux
        cut_e = -11.6
        mask_e = np.array(blazars["Energy_Flux100"]>10**cut_e)
        blazars = blazars[mask_e]

        maps = {
            "RAJ2000": "ra_deg",
            "DEJ2000": "dec_deg",
        }
        
        blazars = blazars.rename(columns=maps)
        blazars.insert(2, 'dec_rad', blazars['dec_deg']*np.pi/180.)
        blazars.insert(2, 'ra_rad', blazars['ra_deg']*np.pi/180.)

        # Apply cut on latitude 
        new_index = np.arange(len(blazars))
        blazars = blazars.set_index(new_index)
        
        mask_GP = [
            is_outside_GP(
                blazars.at[i, 'ra_deg'],blazars.at[i, 'dec_deg']
            ) for i in range(len(blazars))
        ]
        blazars = blazars[mask_GP]
        new_index = np.arange(len(blazars))
        blazars = blazars.set_index(new_index)
        blazars.insert(
            len(blazars.keys()), 'bkg_pdf', np.empty(len(blazars))
        )

        logging.info("Found {0} sources in total".format(len(blazars)))

        return blazars

    @staticmethod
    def set_gp_threshold():
        """Set a cut in galactic latitude for the catalogue
        (exclude the sources with a smaller latitude in absolute value).
        """
        return 10.


class AverageFluxWeightHypothesis(Hypothesis):
    """Hypothesis of constant emission from the fermi blazars
    (average flux as weight)
    """
    name = "average_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        """Weight the catalogue

        Parameters
        ----------
        cat_data: `pandas.DataFrame`
            catalogue to weight
        """
        try:
            return cat_data["Energy_Flux100"]
        except:
            return cat_data['X band map']

        
class BrightestFluxWeightHypothesis(Hypothesis):
    """Hypothesis of constant emission from the fermi blazars
    (average flux as weight). Option to select only the 100 brightest
    [Probably necessary for older tests. Should it be kept?]
    """
    name = "brightest_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        """Weight the catalogue

        Parameters
        ----------
        cat_data: `pandas.DataFrame`
            catalogue to weight
        """
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
        """Weight the catalogue

        Parameters
        ----------
        cat_data: `pandas.DataFrame`
            catalogue to weight
        nu_at: `float`
            neutrino arrival time
        ignore_times: `bool`
            It has no function in here, but maybe keep it
            for compatibility reasons?
        """

        weights = [data_lc[name][nu_at] for name in cat_data['Source_Name']]
                
        return weights
    