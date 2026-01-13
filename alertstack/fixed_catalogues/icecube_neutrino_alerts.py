from alertstack import (
    PointSource, FixedCatalogue, alertstack_data_dir
)
from astropy.io import fits
#import healpy as hp
import logging
import mhealpy as mhp
import numpy as np
import os
import pickle
import resource
from scipy.stats import norm
from scipy import sparse
import time


class NeutrinoAlert(PointSource):

    def __init__(self, ra_deg, dec_deg, time_mjd, weight=1.):
        self.ra_rad = np.radians(float(ra_deg))
        self.dec_rad = np.radians(float(dec_deg))
        self.time_mjd = time_mjd
        self.weight = weight


class CircularisedNeutrinoAlert(NeutrinoAlert):
    '''
    Definition of all the necessary functions to run the analysis with the circularised neutrino alerts
    '''

    def __init__(self, time_mjd, ra, ra_delta, dec, dec_delta, weight=0.5):
        NeutrinoAlert.__init__(self, ra, dec, time_mjd, weight=weight)
        self.sigma = np.radians(
            np.sqrt(0.25 * (ra_delta[0] ** 2 + ra_delta[1] ** 2 + dec_delta[0] ** 2 + dec_delta[1] ** 2)))

    @staticmethod
    def gaussian(delta, sigma):
        return (1. / (2. * np.pi * sigma ** 2.) *
                 np.exp(-0.5 * (delta / sigma) ** 2.))

    @staticmethod
    def bkg_spatial():
        return 1. / (4. * np.pi)

    def eval_spatial_pdf(self, ra, dec):
        delta = self.angular_distance(
            ra, dec, self.ra_rad, self.dec_rad
        )
        return self.gaussian(delta, self.sigma)

    def simulate_position(self):
        sim_ra = self.ra_rad + norm.rvs(scale=self.sigma)
        if sim_ra > 2*np.pi:
            sim_ra -= 2 * np.pi
        elif sim_ra < 0.:
            sim_ra += 2 * np.pi
        sim_dec = self.dec_rad + norm.rvs(scale=self.sigma)
        sim_dec = np.arcsin(np.sin(sim_dec))
        return sim_ra, sim_dec


class HealpixNeutrinoAlert(PointSource):
    '''
    Class of neutrino alerts that will fill the HealpixNeutrinoCatalogue.
    Contains the multiorder probability maps from IceCat-2.
    '''
    def __init__(
        self,
        fits_path,
    ):
        self.fits_path = fits_path
        logging.info("Loading from {0}".format(self.fits_path))
        info, fitsfile = fits.open(fits_path)
        self.header = fitsfile.header
        skymap = fitsfile.data
        
        self.time_mjd = self.header['MJD-OBS']
        self.runid = self.header['RUNID']
        self.eventid = self.header['EVENTID']
        self.nside = self.header['NSIDE']
        self.area_finest_pixel = mhp.nside2pixarea(self.nside)
        self.max_err = np.max([
            self.header["RA_ERR_PLUS_90"],
            self.header["RA_ERR_MINUS_90"],
            self.header["DEC_ERR_PLUS_90"],
            self.header["DEC_ERR_MINUS_90"]
        ])

        self.uniqs = skymap['UNIQ']
        self.probdensity = skymap["PROBDENSITY"]
        self.moc = mhp.HealpixMap(  # Multi-order map
            data=self.probdensity,
            uniq=self.uniqs,
            density=True
        )
        self.pixels = np.arange(len(self.probdensity))
        self.probs = self.probdensity * mhp.nside2pixarea(
            mhp.uniq2nside(self.uniqs)
        )
        self.probs[np.isnan(self.probs)] = 0.  # avoid crash in simulate_position
        self.n_pixels = mhp.nside2npix(self.nside)

        self.ra_deg = self.header["RA"]
        self.dec_deg = self.header["DEC"]
        self.ra_rad = self.ra_deg * np.pi / 180.
        self.dec_rad = self.dec_deg * np.pi / 180.
        self.weight = self.header['P_ASTRO']
        self.energy = self.header['ENERGY']

    # Get the value of the probability from the neutrino likelihood skymap in the given coordinates 
    def signal_pdf(self, ra, dec):
        colat = np.pi / 2. - dec
        long = ra
        #compl_map = np.zeros(hp.nside2npix(self.nside))
        #compl_map[self.pixels] = self.probs
        probdens = self.moc.get_interp_val(colat, long, lonlat=False)
        return probdens * mhp.nside2pixarea(self.nside)
        # return hp.pixelfunc.get_interp_val(compl_map, colat, long, lonlat=False)

    def bkg_spatial_pdf(self):
        return 1./mhp.nside2npix(self.nside)
        # return 1./hp.nside2npix(self.nside)

    # Get coordinates for a point in a healpix grid
    def extract_ra_dec(self, index):
        nside, nestpix = mhp.uniq2nest(self.uniqs[index])
        (colat, ra) = mhp.pix2ang(nside, nestpix, nest=True)
        # (colat, ra) = hp.pix2ang(self.nside, index)
        dec = np.pi / 2. - colat
        # dec = colat
        return ra, dec

    def eval_spatial_pdf(self, ra, dec):
        return self.signal_pdf(ra, dec) # /self.bkg_spatial_pdf()

    # Simulate a random position (return coordinates) weighted by the neutrino PSF
    def simulate_position(self):
        ind = np.random.choice(a=self.pixels, p=self.probs)
        pos = self.extract_ra_dec(ind)
        return pos

class CircularisedNeutrinoAlertCatalogue(FixedCatalogue):
    '''
        This class contains a subset of alerts that were published for the TXS paper. The catalog 
        only includes circularised errors and signalness = 0.5 for every alert. 
    '''
    @staticmethod
    def parse_data():
        nu_objs = []
        with open(os.path.join(alertstack_data_dir, "catalog_of_alerts.txt"), "r") as f:
            for line in f.readlines():
                if line[0] not in ["#", "\n"]:
                    if "retracted" not in line:
                        vals = [x for x in line.split(" ") if x not in [""]]
                        time = vals[0]
                        ra = vals[1]
                        dec = vals[3]
                        ra_delta = [float(x)/2.5 for x in vals[2][1:-1].split(",")]
                        dec_delta = [float(x)/2.5 for x in vals[4][1:-2].split(",")]
                        ra_delta *= abs(np.cos(dec_delta))

                        nu_objs.append(CircularisedNeutrinoAlert(time, ra, ra_delta, dec, dec_delta))

        return nu_objs


    def add_sim_alerts(self, n):

        nu_objs = []

        for _ in range(n):
            sigma = 0.2
            ra_delta = [sigma, sigma]
            dec_delta = [sigma, sigma]
            time = 0.0
            ra = np.degrees(np.random.uniform() * 2 * np.pi)
            dec = np.degrees(np.arccos(2.*np.random.uniform() - 1) - np.pi/2.)
            nu_objs.append(CircularisedNeutrinoAlert(time, ra, np.array(ra_delta), dec, np.array(dec_delta)))

        self.data += nu_objs

try:
    skymap_dir = os.environ['NU_SKYMAP_DIR']
except KeyError:
    logging.warning("No NU_SKYMAP_DIR variable set. If you do not set this, importing a "
                   "HealpixNeutrinoAlertCatalogue will raise an error.")

class HealpixNeutrinoAlertCatalogue(FixedCatalogue):
    '''
        Catalog containing all the neutrino alerts in the alert catalog v2. It loads the healpix skymaps 
        with the likelihood information from the Millipede scans. 
    '''
    @staticmethod
    def parse_data():
        nu_objs = []

        logging.info("Loading from {0}".format(skymap_dir))
        
        files = [x for x in os.listdir(skymap_dir) if ".multiorder.fits.gz" in x]
        for filename in files:
            #print(filename)
            path =  os.path.join(skymap_dir, filename)
            nu = HealpixNeutrinoAlert(path)
            nu_objs.append(nu)
            
        return nu_objs
