from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
import copy
import healpy as hp
import numpy as np
import os
import logging
from pathlib import Path
import pandas as pd
import pickle as pkl
from alertstack import AnisotropicExtragalacticCatalogue, Hypothesis, is_outside_GP, alertstack_data_dir
from numpy.lib.recfunctions import rename_fields
from astropy import units as u
import astropy.io.ascii


class FlairesCatalogue(AnisotropicExtragalacticCatalogue):
    '''
    Loads the 528 accretion flares used in Flairestack.
    '''

    @staticmethod
    def parse_data():

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        # Load catalog from paper J. Necker et al. (2025)
        names = [
            "Source_Name", "ra_deg", "dec_deg", "HW1mag", "HW2mag",
            "NEWS", "AllWISE", "PS1", "AllWISEcntr", "WISE",
            "AllWISEDes", "NEDLVSIndex", "NEDLVSName", "ParentSamp", "z",
            "e_z", "r_z", "SDSSdist", "SDSSclass", "TNSdist",
            "TNSobjtype1", "TNSobjtype2",
            "TNSname", "TNSdate", "milliqdist", "milliqtype", "MirongName",
            "WTPName", "RefTime", "x2W1", "npointsW1", "FmedW1",
            "x2W2", "npointsW2", "FmedW2", "FbslW1", "e_FbslW1",
            "FbslW2", "e_FbslW2", "startW1", "endW1", "endedW1",
            "startW2", "endW2", "endedW2", "strengthW1", "strengthW2",
            "varW1", "varW2", "MaxFluxW1", "e_MaxFluxW1", "MaxFluxW2",
            "e_MaxFluxW2", "FluenceW1", "e_FluenceW1", "FluenceW2", "e_FluenceW2",
            "Sep", "PeakLbol", "PeakTime", "Ebol", "Fluencebol",
        ]
        
        flaires_df = pd.read_fwf(
            os.path.join(alertstack_data_dir, "flaires.dat"),
            names = names,
        )
        flaires_df.insert(1, 'dec_rad', flaires_df['dec_deg']*np.pi/180.)
        flaires_df.insert(1, 'ra_rad', flaires_df['ra_deg']*np.pi/180.)
        flaires_df.insert(
            5,
            't-peak',
            flaires_df["PeakTime"].to_numpy()*(1+flaires_df["z"].to_numpy())+flaires_df["RefTime"].to_numpy()
        )
        flaires_df.insert(len(names)+3, 'bkg_pdf', np.empty(len(flaires_df)))

        # Select the same sources as in Flairestack.
        # Exclude all sources without a bolometric luminosity.
        fluencemask = ~np.isnan(flaires_df['Fluencebol'])
        # Exclude all sources identified as synchrotron emitters from jets.
        syncrojetmask = ~(
            (flaires_df["milliqtype"]=="BR")  |
            (flaires_df["milliqtype"]=="BRX") |
            (flaires_df["milliqtype"]=="KRX") |
            (flaires_df["milliqtype"]=="QR")  |
            (flaires_df["milliqtype"]=="QRX") |
            (flaires_df["milliqtype"]=="Q")   |
            (flaires_df["milliqtype"]=="QX")  |
            (flaires_df["milliqtype"]=="qRX") |
            (flaires_df["milliqtype"]=="q")   | 
            (flaires_df["milliqtype"]=="qR")
        )

        # Cut on declination and GP, to be removed in future
        decmask = flaires_df["dec_deg"] >= -30.
        eq_coords = SkyCoord(
            ra=np.array(flaires_df["ra_rad"])*u.rad,
            dec=np.array(flaires_df["dec_rad"])*u.rad,
            frame='icrs',
        )
        gal_lats = eq_coords.galactic.b.deg
        gallatmask = np.abs(gal_lats) > 8.

        
        timemask = flaires_df['t-peak'] > 55695  # Day of the first alert
        
        acflaresmask = np.array(
            fluencemask & syncrojetmask & timemask  # & decmask & gallatmask
        )
        acflares_df = flaires_df[acflaresmask]

        new_index = np.arange(len(acflares_df))
        acflares_df = acflares_df.set_index(new_index)

        print("Found {0} sources in total".format(len(acflares_df)))
        logging.info("Found {0} sources in total".format(len(acflares_df)))

        return acflares_df

    @staticmethod
    def set_gp_threshold():
        return 0.

    @staticmethod
    def set_min_declination():
        # Minimal declination is set to the minimum of ZTF.
        return -90.

    @staticmethod
    def set_nside():
        # nside of the bkg distribution.
        return 128

    def set_bkg_distribution(self):
        return  self.generate_bkg_distribution_allsky(
                self.data, nside=16, hd_nside=128, sigma_smoothing=8.,
        )

class FluencebolHypothesis(Hypothesis):
    name = "bolometric_fluence_weight"

    @staticmethod
    def weight_catalogue(cat_data, nu_at, ignore_times=False):
        """
        Consider the bolometric luminosity [mJy s-1] as weight.
        The time window consits of 1 year before the IR peak.
        If the flare is not within 1 year after the neutrino the
        weight is set to zero.
        """
        # This option is for the injections where we do not care about
        # selecting with the time window
        if ignore_times:
            return cat_data["Fluencebol"]
        
        # Consider the coincidence only if within 1 year from neutrino
        weights = np.zeros(len(cat_data))
        peak_times = cat_data["t-peak"]
        weightmask = (
            (nu_at - peak_times <= 0.) &
            (nu_at - peak_times >= -365.)
        )
        weights[weightmask] = cat_data["Fluencebol"][weightmask]
        return weights
    

