import numpy as np
import logging
import pickle
from scipy.stats import norm
from scipy import sparse
from alertstack import PointSource, FixedCatalogue, is_outside_GP
import healpy as hp
import os


class NeutrinoAlert(PointSource):

    def __init__(self, ra_deg, dec_deg, time_mjd, weight=1.):
        self.ra_rad = np.radians(float(ra_deg))
        self.dec_rad = np.radians(float(dec_deg))
        self.time_mjd = time_mjd
        self.weight = weight


class CircularisedNeutrinoAlert(NeutrinoAlert):

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

    def __init__(self, pkl_path):
        self.pkl_path = pkl_path
        logging.info("Loading from {0}".format(pkl_path))
        with open(self.pkl_path, "rb") as f:
            self.pkl_dict = pickle.load(f)
        self.mask = sparse.load_npz(self.pkl_dict["output_path"]).toarray()[0]
        self.probs = np.zeros(len(self.mask))
        self.probs[self.mask] = np.array(self.pkl_dict["prob"])
        self.n_pixels = float(len(self.probs))
        self.nside = hp.pixelfunc.npix2nside(self.n_pixels)
        self.ra_rad,self.dec_rad = self.extract_ra_dec(np.where(self.probs == np.max(self.probs)))
        self.ra_deg = self.pkl_dict["RA"]  
        self.dec_deg = self.pkl_dict["DEC"]

        try:
            self.weight = self.pkl_dict["SIGNAL"]
            if type(self.weight)== str: # some neutrinos don't have signalness? just 2 or 3
                self.weight = 0.5
        except KeyError:
            self.weight = 0.5

    def signal_pdf(self, ra, dec):
        colat = np.pi / 2. - dec
        return hp.pixelfunc.get_interp_val(self.probs, colat, ra, lonlat=False)

    def bkg_spatial_pdf(self):
        return 1./self.n_pixels

    def extract_ra_dec(self, index):
        (colat, ra) = hp.pix2ang(self.nside, index)
        dec = np.pi / 2. - colat
        return ra, dec

    def eval_spatial_pdf(self, ra, dec):
        return self.signal_pdf(ra, dec)/self.bkg_spatial_pdf()

    def simulate_position(self):
        ind = np.random.choice(int(self.n_pixels), p=self.probs)
        pos = self.extract_ra_dec(ind)
        return pos

class CircularisedNeutrinoAlertCatalogue(FixedCatalogue):

    @staticmethod
    def parse_data():
        nu_objs = []
        with open("data/catalog_of_alerts.txt", "r") as f:
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

#dir = "/Users/robertstein/Realtime_Stuff/alert_archive/output_raw_fits/compressed_files/"
dir = "/Users/crislagual/Documents/phd/1_year/alertstack/data/alerts_archive/compressed_files"

class HealpixNeutrinoAlertCatalogue(FixedCatalogue):

    @staticmethod
    def parse_data():
        nu_objs = []

        logging.info("Loading from {0}".format(dir))

        files = [x for x in os.listdir(dir) if ".pkl" in x]
        for filename in files:
            path =  os.path.join(dir, filename)
            nu = HealpixNeutrinoAlert(path)
            nu_objs.append(nu)
            #    if is_outside_GP(np.rad2deg(nu.ra_rad),np.rad2deg(nu.dec_rad)): # only store those outside GP
            #        nu_objs.append(nu)

        return nu_objs

HealpixNeutrinoAlertCatalogue()

