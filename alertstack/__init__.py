import numpy as np
import healpy as hp
import random
from astropy import units as u
from astropy.coordinates import SkyCoord
import os
import pickle
import copy

import time

alertstack_data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data/")

cat_dtype = np.dtype([
    ("Ra", float),
    ("Dec", float)
])


class PointSource:

    def eval_spatial_pdf(self, ra, dec):
        raise NotImplementedError

    @staticmethod
    def angular_distance(lon1, lat1, lon2, lat2):
        """calculate the angular distince along the great circle
        on the surface of a sphere between the points
        (`lon1`,`lat1`) and (`lon2`,`lat2`)
        This function Works for equatorial coordinates
        with right ascension as longitude and declination
        as latitude. This function uses the Vincenty formula
        for calculating the distance.
        Parameters
        ----------
        lon1 : array_like
          longitude of first point in radians
        lat1 : array_like
          latitude of the first point in radians
        lon2 : array_like
          longitude of second point in radians
        lat2 : array_like
          latitude of the second point in radians
        """
        c1 = np.cos(lat1)
        c2 = np.cos(lat2)
        s1 = np.sin(lat1)
        s2 = np.sin(lat2)
        sd = np.sin(lon2 - lon1)
        cd = np.cos(lon2 - lon1)

        return np.arctan2(
            np.hypot(c2 * sd, c1 * s2 - s1 * c2 * cd),
            s1 * s2 + c1 * c2 * cd
        )

    def simulate_position(self):
        raise NotImplementedError

    def eval_source_weight(self):
        return self.weight

class Catalogue:

    def __init__(self):
        self.data = self.parse_data()
        self.gp_threshold = self.set_gp_threshold()
        self.min_declination = self.set_min_declination()
        self.nside = self.set_nside()
        self.npix = self.set_npix()
        self.numap_nside = 1024
        self.numap_npix = hp.nside2npix(self.numap_nside)
        self.bkg_distribution = self.set_bkg_distribution()
        self.set_bkg_pdf_per_source(self.data)

    @staticmethod
    def parse_data():
        return NotImplementedError

    @staticmethod
    def set_gp_threshold():
        return NotImplementedError

    @staticmethod
    def set_min_declination():
        return NotImplementedError

    @staticmethod
    def set_nside():
        return NotImplementedError

    @staticmethod
    def generate_bkg_distribution_allsky(
        cat, nside, hd_nside=128, sigma_smoothing=15.,
    ):
        """
        Given the coordinates of the sources, it returns an all-sky
        probability map.
    
        :param catalog: catalog with info regarding coordinates
        :param nside: nside for healpix histogramming of sources
        :param hd_nside: up to which nside the map must be upgraded
        :param sigma_smoothing [deg]: smoothing to apply to the 
        healpix histogram
        :return: the healpix map with the probabilities
        """

        theta = np.pi/2. - cat["dec_rad"].to_numpy()
        phi = cat["ra_rad"].to_numpy()
        
        bins_per_source = hp.ang2pix(nside, theta, phi)
        bins = np.arange(hp.nside2npix(nside)+1)
        counts_per_bin, _ = np.histogram(bins_per_source, bins=bins)
        bins_probs = counts_per_bin/np.sum(counts_per_bin)
        bins_probs = hp.ud_grade(bins_probs, hd_nside)
        bins_probs = bins_probs / np.sum(bins_probs)
        bins_probs = hp.sphtfunc.smoothing(
            bins_probs, sigma=sigma_smoothing*np.pi/180.
        )
        bins_probs[bins_probs<0.] = 0.
        bins_probs = bins_probs / np.sum(bins_probs)
        return bins_probs

    @staticmethod
    def apply_cuts_on_bkg_distribution(
        bins_probs, min_declination=-25, gp_threshold=8.
    ):
        """
        Upgrades the map resolution and applies the necessary cuts.
    
        :param bins_probs: the initial healpix map (an array)
        :param min_theta [deg]: lower cut on declination
        :param max_gal_lat [deg]: exclude the galactic plane up to this
        latitude
        :return: upgraded map with cuts applied
        """

        hd_npix = len(bins_probs)
        hd_nside = hp.npix2nside(hd_npix)
        hd_cothetas, hd_phis = hp.pix2ang(
            hd_nside, np.arange(hd_npix)
        )
        hd_thetas = np.pi/2. - hd_cothetas
        low_theta_mask = hd_thetas < min_declination * np.pi / 180.
        hd_icrs = SkyCoord(
            ra=hd_phis*u.rad,
            dec=hd_thetas*u.rad,
            frame='icrs',
            unit='rad',
        )
        hd_bs = hd_icrs.galactic.b.deg
        galactic_mask = np.abs(hd_bs) < gp_threshold
        
        bins_probs[low_theta_mask | galactic_mask] = 0.
        bins_probs = bins_probs / np.sum(bins_probs)
    
        return bins_probs

    def set_npix(self):
        return NotImplementedError

    def bkg_spatial_pdf(self, ra, dec):
        # Evaluated the background probability in a specified direction
        # and returns the value that corresponds to a map with the same
        # nside as a neutrino map.
        bkg_pix = hp.ang2pix(self.nside, np.pi/2. - dec, ra)
        numap_npix = hp.nside2npix(self.numap_nside)
        return self.bkg_distribution[bkg_pix] * self.npix / numap_npix

    def set_bkg_distribution(self):
        return NotImplementedError

    def set_bkg_pdf_per_source(self, cat):
        return NotImplementedError

class FixedCatalogue(Catalogue):

    def __getitem__(self, item):
        return self.data[item]

    def __iter__(self):
        return self.data.__iter__()



class ScrambleCatalogue(Catalogue):

    def unblind(self):
        return self.data

    def scramble(self):
        return NotImplementedError

    @staticmethod
    def extract_ra_dec(nside, index):
        (colat, ra) = hp.pix2ang(nside, index, nest=True)
        dec = np.pi / 2. - colat
        return ra, dec

    def return_ra_dec(self):
        return NotImplementedError

class IsotropicExtragalacticCatalogue(ScrambleCatalogue):

    def __init__(self, nside=1024):
        ScrambleCatalogue.__init__(self)

    def scramble_positions_outside_GP(self, gp_cut=10., min_dec_deg=-90.):
        # Scramble positions directly outside the galactic plane

        def perform_scramble_outside_GP(data=self.data):
            
            gal_l_vals = np.random.uniform(low=0, high=2*np.pi, size=len(data))
         
            # divide the sky into two areas, half of the blazars on each side
            threshold = np.deg2rad(gp_cut) # 10 degrees 
            size_half = int(len(data) / 2)
    
            gal_b_vals_up = np.arccos(
                2*np.random.uniform(
                    low=0,high=0.5*(1-np.sin(threshold)),size=size_half
                )-1
            ) - np.pi/2.
            gal_b_vals_down = np.arccos(
                2 * np.random.uniform(
                    low=0.5*(1+np.sin(threshold)),
                    high=1,
                    size=(len(data) - size_half)
                ) - 1
            ) - np.pi / 2.
            gal_b_vals = np.concatenate(
                (gal_b_vals_down,gal_b_vals_up),axis=0
            )
            np.random.shuffle(gal_b_vals)
    
            gal = SkyCoord(
                l = gal_l_vals*u.rad, b = gal_b_vals*u.rad, frame='galactic'
            )
            ra_vals = gal.icrs.ra.rad
            dec_vals = gal.icrs.dec.rad

            return ra_vals, dec_vals

        ra_vals, dec_vals = perform_scramble_outside_GP()
        too_low_mask = dec_vals * 180. / np.pi < min_dec_deg

        while np.sum(too_low_mask) > 0:
            # repeat scramble only for sources too low in declination
            (
                ra_vals[too_low_mask], dec_vals[too_low_mask]
            ) = perform_scramble_outside_GP(self.data[too_low_mask])
            too_low_mask = dec_vals * 180. / np.pi < min_dec_deg
        

        return ra_vals, dec_vals
    
    def scramble_positions(self, min_dec=-90.):
        # Scramble positions
        
        ra_vals = np.random.uniform(size=len(self.data)) * 2 * np.pi
        min_dec_rad = min_dec * np.pi / 180.
        dec_vals = np.arccos(
            (
                1-np.sin(min_dec_rad)
            )*np.random.uniform(size=len(self.data)) + np.sin(min_dec_rad)
        ) - np.pi/2.
        return ra_vals, dec_vals


def is_outside_GP(ra,dec, threshold=10.0):
    # Check if a position in the sky (in degrees) is outside of the galactic plane

    eq = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    gal = eq.galactic
    threshold_GP = threshold*u.deg

    return abs(gal.b) > threshold_GP

class Hypothesis:
    name = None

    def __init__(self, fixed_catalogue, min_E=0.):
    #  min_E added to test minimum sensitive energy
        self.fixed_catalogue = fixed_catalogue

        nu_energies = np.array([nu.energy for nu in fixed_catalogue])
        energymask = nu_energies >= min_E
        
        
        if self.name == 'strength_flux_weight':
            # Select only neutrinos which can be coincident with the 63 accretion flares
            nutimes = np.array(
                [nu.time_mjd for nu in fixed_catalogue]
            )
            nudecs = np.array(
                [nu.dec_deg for nu in fixed_catalogue]
            )
            nudecsplus = np.array(
                [nu.header["DEC_ERR_PLUS_90"] for nu in fixed_catalogue]
            )
            nutimes = np.array(nutimes)
            minflarestime = 58261.4
            maxflarestime = 58978.3 + 365.
            timemask = (
                (nutimes >= minflarestime) & (nutimes <= maxflarestime)
            )
            spatialmask = (nudecs + nudecsplus) > -25.
            self.fixed_catalogue = np.array(fixed_catalogue.data)[timemask & spatialmask & energymask]
            #selected_nus = [f"{nu.header["RUNID"]} {nu.header["EVENTID"]}" for nu in self.fixed_catalogue]
            #selected_nus.sort()
            #for i, nu in enumerate(selected_nus):
            #    print(i, nu)
        elif self.name == 'bolometric_fluence_weight':
            # Select only neutrinos that can be coincident with the 524 accretion flares
            nutimes = np.array(
                [nu.time_mjd for nu in fixed_catalogue]
            )
            nutimes = np.array(nutimes)
            minflarestime = 55465.79 - 365.
            maxflarestime = 59747.02
            timemask = (
                (nutimes >= minflarestime) & (nutimes <= maxflarestime)
            )
            self.fixed_catalogue = np.array(fixed_catalogue.data)[timemask & energymask]
            #selected_nus = [f"{nu.header["RUNID"]} {nu.header["EVENTID"]}" for nu in self.fixed_catalogue]
            #selected_nus.sort()
            #for i, nu in enumerate(selected_nus):
            #    print(i, nu)
        else:
            self.fixed_catalogue = np.array(fixed_catalogue.data)[energymask]
            
        self.source_weights = np.array([source.eval_source_weight() for source in self.fixed_catalogue])

    @staticmethod
    def weight_catalogue(cat_data):
        return NotImplementedError

    
    def calculate_llh(self, cat_data, savedata=None, gp_threshold=10.0):
        '''
        Calculate the TS_i of each neutrino as TS_i = log(S/B), where S = max(S_spatial * signalness * w_blazar) 
        and B = B_spatial. If the neutrino is in the Galactic Plane or TS_i < 0, then TS_i = 0 (S/B = 1, 
        choose background hypothesis). The final TS of the trial is simply TS = sum(TS_i).
        '''
            
        cat_mask  = is_outside_GP(
            np.array(np.rad2deg(cat_data["ra_rad"])),
            np.array(np.rad2deg(cat_data["dec_rad"])),
            threshold = gp_threshold,
        ) == False

        if savedata is not None:
            final = []

        if (
            (self.name != 'monthly_flux_weight') and 
            (self.name != 'strength_flux_weight') and
            (self.name != 'bolometric_fluence_weight')
        ):
            cat_weights = self.weight_catalogue(cat_data) # w_blazars
            density = np.sum(cat_weights)

        lh_array = 0.

        for i, source in enumerate(self.fixed_catalogue): # loop over neutrinos
            
            max_dist = 4 * source.max_err * np.pi/180.
            
            s_ra, s_de = source.ra_rad, source.dec_rad
            corads = copy.copy(cat_data["ra_rad"]) - s_ra + np.pi
            corads %= 2*np.pi

            dist_mask = np.logical_and(
                np.abs(s_de - cat_data["dec_rad"]) < max_dist,
                np.abs(corads - np.pi) < max_dist
            )
            
            spatial_pdf = np.zeros(len(cat_data["dec_rad"]))

            spatial_pdf[dist_mask] = source.eval_spatial_pdf(
                cat_data["ra_rad"][dist_mask],
                cat_data["dec_rad"][dist_mask]
            ) / cat_data["bkg_pdf"][dist_mask] # * (4 * np.pi)


            if self.name != 'average_radio_flux_weight':
                spatial_pdf_mask = np.where(
                    cat_mask,
                    0.0,
                    np.array(spatial_pdf)
                ) # if the neutrino is in the GP, TS_i = 0
            else:
                spatial_pdf_mask = spatial_pdf # don't need to mask the GP for the radio catalog


            source_weight = self.source_weights[i] # signalness


            if (
                (self.name == 'monthly_flux_weight') or 
                (self.name == 'strength_flux_weight') or
                (self.name == 'bolometric_fluence_weight')
            ):
                cat_weights = self.weight_catalogue(cat_data, source.time_mjd) # w_blazars
                density = np.sum(cat_weights)

            if density == 0.:
                prob = 1.
            else:
                prob = max(source_weight * spatial_pdf_mask * cat_weights / density) # S/B
            if prob < 1.:
                prob = 1.

            lh_array += np.log(prob) # TS = log(S/B)
            
            if savedata is not None:
                ind = np.argmax(source_weight * spatial_pdf_mask * cat_weights / density)
                final.append([source.pkl_path, cat_data[ind]['Source_Name'], np.log(prob)])
            
        if savedata is not None:
            with open(os.path.join(savedata,"correlations.pkl"), "wb") as fp:
                pickle.dump(final, fp)

        llh = lh_array
        return llh


    def inject_signal(self, cat, fraction):
        '''
        Create signal trials by injecting correlations between neutrino alerts and the catalog sources.
        '''


        nucat = self.fixed_catalogue

        
        n_exp = fraction * np.sum(np.array(self.source_weights)) # Choose expected number of neutrinos (astrophysical or not) to have correlations
        n_inj = np.random.poisson(n_exp) # Get number of neutrinos with correlations (poisson fluctuation)

        if n_inj > len(cat):
            raise Exception("Trying to inject more sources than there are entries in the catalogue! \n"
                            "There are {0} entries in the catalogue, and the expectation for injection is {1}. \n"
                            "`Applying random poisson noise, we are trying to inject {2} this trial".format(
                len(cat), n_exp, n_inj
            ))

        if n_inj > 0:

            # Choose which neutrinos will have a counterpart (each neutrino can only be injected once in each trial)
            ind = np.random.choice(
                len(np.array(self.source_weights)),
                size=n_inj, 
                p=np.array(self.source_weights)/np.sum(np.array(self.source_weights)),
                replace=False
            )
            
            inj_cat = []

            inj_sources = []

            for i in ind:
                fixed_source = np.array(self.fixed_catalogue.data)[i]

                # Choose which counterpart, according to the weighting scheme
                if (
                    (self.name == 'monthly_flux_weight') or 
                    (self.name == 'strength_flux_weight') or
                    (self.name == 'bolometric_fluence_weight')
                ):
                    weights = self.weight_catalogue(
                        cat,
                        fixed_source.time_mjd,
                        ignore_times=True,
                    )
                else:
                    weights = self.weight_catalogue(cat)
                weights /= np.sum(weights)

                source_names = cat['Source_Name'].to_numpy()

                # Mask already chosen sources
                j = np.random.choice(len(weights), p=weights)
                while source_names[j] in inj_sources:
                    j = np.random.choice(len(weights), p=weights)
                        
                inj_sources.append(source_names[j])

                # Simulate new source position
                fixed_source.probs = fixed_source.probs/ np.sum(fixed_source.probs)
                cat.at[j, 'ra_rad'], cat.at[j, 'dec_rad'] = fixed_source.simulate_position()

                # Change the date so that it is coincident with the neutrino
                if (
                    (self.name == 'strength_flux_weight')
                ):
                    new_time = fixed_source.time_mjd - np.random.random()*365.
                    cat.at[j, 't-peak'] = new_time
                elif (
                    (self.name == 'bolometric_fluence_weight')
                ):
                    new_time = fixed_source.time_mjd + np.random.random()*365.
                    cat.at[j, 't-peak'] = new_time
                    

        return cat


class UniformPriorHypothesis(Hypothesis):
    name = "uniform_prior"

    @staticmethod
    def weight_catalogue(cat_data):
        return np.ones(len(cat_data))
