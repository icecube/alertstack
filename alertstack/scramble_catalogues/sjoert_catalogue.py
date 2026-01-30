from astropy.io import fits
import copy
import healpy as hp
import numpy as np
import os
import logging
from pathlib import Path
import pandas as pd
import pickle as pkl
from alertstack import AnisotropicExtragalacticCatalogue, Hypothesis, is_outside_GP, alertstack_data_dir
from hellolancel.TS_input import (
    p_flux_bg, p_flux_sig, p_strength_bg, p_strength_sig
)
from numpy.lib.recfunctions import rename_fields
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.io.ascii


class AccretionFlaresSjoertCatalogue(AnisotropicExtragalacticCatalogue):
    '''Loads Sjoert's 63 accretion flares.
    '''

    @staticmethod
    def parse_data():
        """Load the catalogue.
        """

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        # Load catalog from paper S. van Velzen et al. (2024)
        # Here without R.A. and Dec. To add from other data released
        # on Zenodo (https://zenodo.org/records/7026636)
        acflares_df = pd.read_csv(
            os.path.join(alertstack_data_dir, "sjoert_catalog.txt"),
            sep='\t',
            header=0,
            names=[
                'Source_Name', 't-peak', 'Riseb', 'Fadeb', 'ΔFIR/Frms',
                'ΔFIR', 'z', 'MBH', 'PAGN', 'Spectro.',
            ]
        )

        # Data release on Zenodo with additional information
        # regarding R.A. and Dec..
        ztf_acflares = astropy.io.ascii.read(
            os.path.join(
                alertstack_data_dir, "ZTF_neoWISE_flares_acflares.dat"
            )
        )

        # Since the same names of the paper are not reported,
        # recognize the correct flares through MJD.
        mjds = np.array(ztf_acflares['flare_peak_jd']) - 2400000.5
        mjds = np.array([round(d,1) for d in mjds])
        ras = []
        des = []
        mjds_tdes = acflares_df['t-peak'].to_numpy()
        for index, mjd in enumerate(mjds_tdes):
            if mjd in mjds:
                mjd_index = np.where(mjds==mjd)[0][0]
                ra = ztf_acflares['ra'][mjd_index]
                de = ztf_acflares['dec'][mjd_index]
                ras.append(ra)
                des.append(de)
        ras = np.array(ras)
        des = np.array(des)
        acflares_df.insert(1, 'dec_deg', des)
        acflares_df.insert(1, 'ra_deg', ras)
        acflares_df.insert(1, 'dec_rad', des*np.pi/180.)
        acflares_df.insert(1, 'ra_rad', ras*np.pi/180.)

        # Determine once weights to speed up the code
        logstrength = np.log10(np.array(acflares_df["ΔFIR/Frms"]))
        echo_fluxes = np.array(
            [float(flux.split("±")[0])*1e-3 for flux in acflares_df['ΔFIR']]
        )
        logflux = np.log10(echo_fluxes)
        weightstrength = p_strength_sig(logstrength) / p_strength_bg(logstrength)
        weightflux = p_flux_sig(logflux) / p_flux_bg(logflux)
        acflares_df.insert(6, 'weight', weightstrength * weightflux)
        acflares_df.insert(7, 'bkg_pdf', np.empty(len(acflares_df)))

        logging.info("Found {0} sources in total".format(len(acflares_df)))

        return acflares_df

    @staticmethod
    def set_gp_threshold():
        """Cut of 8 deg same as in paper S. van Velzen at al. (2024).
        """
        return 8.

    @staticmethod
    def set_min_declination():
        """Minimal declination is set to the minimum of ZTF.
        """
        return -25.

    @staticmethod
    def set_nside():
        """nside of the bkg distribution.
        """
        return 128

    def set_bkg_distribution(self):
        """Contains the logic necessary to generate an appropriate
        background distribution for the specific catalogue.
        """
        return self.apply_cuts_on_bkg_distribution(
            self.generate_bkg_distribution_allsky(
                self.data, nside=4, hd_nside=self.nside,
            ),
            self.min_declination,
            self.gp_threshold
        )


class StrengthFluxWeightHypothesis(Hypothesis):
    """Class for the hypothesis of neutrino emission proportional
    to the probability ratio used in the paper S. van Velzen et al. (2024).
    """
    name = "strength_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data, nu_at, ignore_times=False):
        """
        Get same exact weights for the sources as the ones used in the
        paper S. van Velzen et al. (2024)

        Parameters
        ----------
        cat_data: `pandas.DataFrame`
            catalogue to weight
        nu_at: `float`
            Neutrino arrival time, in the form of modified julian date.
        ignore_times: `bool`
            If True, ignore the time window and weight only with the
            bolometric fluence.
        """

        # This option is for the injections where we do not care about
        # selecting with the time window
        if ignore_times:
            return cat_data["weight"]
            
        # Consider the coincidence only if within 1 year from ZTF peak
        cat_times = cat_data["t-peak"].to_numpy()
        weights = np.zeros(len(cat_data))
        weightmask = (
            (nu_at - cat_times <= 365.) &
            (nu_at - cat_times >= 0.)
        )
        weights[weightmask] = cat_data["weight"][weightmask]

        return weights
    

