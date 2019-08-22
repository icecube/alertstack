import numpy as np
from scipy.stats import norm
from alertstack import PointSource, FixedCatalogue

class NeutrinoAlert(PointSource):

    def __init__(self, ra_deg, dec_deg, time_mjd, weight=1.):
        self.ra_rad = np.radians(float(ra_deg))
        self.dec_rad = np.radians(float(dec_deg))
        self.time_mjd = time_mjd
        self.weight = weight

    def eval_source_weight(self):
        return self.weight

class CircularisedNeutrinoAlert(NeutrinoAlert):

    def __init__(self, time, ra, ra_delta, dec, dec_delta, weight=0.5):
        NeutrinoAlert.__init__(self, ra, dec, time, weight=weight)
        self.sigma = np.sqrt(0.25 * (ra_delta[0] ** 2 + ra_delta[1] ** 2 + dec_delta[0] ** 2 + dec_delta[1] ** 2))

    @staticmethod
    def gaussian(delta, sigma):
        return norm.pdf(x=delta, scale=sigma) ** 2

    @staticmethod
    def bkg_spatial():
        return 1. / (4. * np.pi)

    def eval_spatial_pdf(self, ra, dec):
        delta = self.angular_distance(
            ra, dec, self.ra_rad, self.dec_rad
        )
        return self.gaussian(delta, self.sigma) / self.bkg_spatial()


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
                        ra_delta = [float(x) for x in vals[2][1:-1].split(",")]
                        dec_delta = [float(x) for x in vals[4][1:-2].split(",")]
                        nu_objs.append(CircularisedNeutrinoAlert(time, ra, ra_delta, dec, dec_delta))

        return nu_objs

    def __iter__(self):
        return self.data.__iter__()