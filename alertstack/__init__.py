import numpy as np
import healpy as hp
from scipy.stats import norm
import random


cat_dtype = np.dtype([
    ("Ra", np.float),
    ("Dec", np.float)
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

class Catalogue:

    def __init__(self):
        self.data = self.parse_data()

    @staticmethod
    def parse_data():
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


gal_plane_1024 = []


class IsotropicExtragalacticCatalogue(ScrambleCatalogue):

    def __init__(self, nside=1024):
        ScrambleCatalogue.__init__(self)
        # self.nside = nside
        # print("NSIDE = {0}, Max Pixel Radius = {1} deg".format(nside, np.degrees(hp.max_pixrad(nside))))
        # self.cone_ids = [x for x in range(hp.nside2npix(self.nside)) if x not in gal_plane_1024]

    def scramble_positions(self):

        ra_vals = np.random.uniform(size=len(self.data)) * 2 * np.pi
        dec_vals = np.arccos(2.*np.random.uniform(size=len(self.data)) - 1) - np.pi/2.
        return ra_vals, dec_vals
        # indexes = np.random.choice(len(self.cone_ids), size=len(self.data))
        # return self.extract_ra_dec(self.nside, indexes)

class Hypothesis:
    name = None

    def __init__(self, fixed_catalogue):
        self.fixed_catalogue = fixed_catalogue
        self.source_weights = np.array([source.eval_source_weight() for source in self.fixed_catalogue])
        # self.source_weights /= np.mean(self.source_weights)

    @staticmethod
    def weight_catalogue(cat_data):
        return NotImplementedError

    def calculate_llh(self, cat_data):
        cat_weights = self.weight_catalogue(cat_data)
        density = np.sum(cat_weights)# / (4 * np.pi)
        # lh_array = np.zeros(len(cat_data))
        lh_array = 0.
        for i, source in enumerate(self.fixed_catalogue):
            spatial_pdf = source.eval_spatial_pdf(cat_data["ra_rad"], cat_data["dec_rad"]) * (4 * np.pi)
            source_weight = self.source_weights[i]

            prob = np.sum(source_weight * spatial_pdf * cat_weights / density)

            # print(prob)
            # input("?")
            #
            # if prob > 0:
            lh_array += np.log(prob + 1.)

        llh = lh_array - np.log(np.sum(self.source_weights))
        # llh = np.log(np.sum(lh_array) + 1.) - np.log(np.sum(self.source_weights))# - np.log(np.sum(cat_weights))
        return llh

    # def calculate_likelihood(self, cat_data, source):
    #     cat_weights = self.weight_catalogue(cat_data)
    #     spatial_weights = source.eval_spatial_pdf(cat_data)
    #     return np.sum(cat_weights)
    #
    # def stack_llh(self, hypo, cat):
    #     lh = 0
    #     for source in self.fixed_catalogue:
    #         lh += self.calculate_likelihood(hypo, cat, source)
    #     return np.log(lh)

    def inject_signal(self, cat, fraction):

        n_exp = fraction * float(len(self.source_weights))
        n_inj = np.random.poisson(n_exp)

        if n_inj > len(cat):
            raise Exception("Trying to inject more sources than there are entries in the catalogue! \n"
                            "There are {0} entries in the catalogue, and the expectation for injection is {1}. \n"
                            "`Applying random poisson noise, we are trying to inject {2} this trial".format(
                len(cat), n_exp, n_inj
            ))

        if n_inj > 0:

            # Choose which fixed source will have a counterpart
            ind = np.random.choice(len(self.source_weights), size=n_inj,
                                  p=self.source_weights/np.sum(self.source_weights))

            inj_cat = []

            inj_sources = []

            for i in ind:
                fixed_source = self.fixed_catalogue[i]

                mask = np.array([k not in inj_sources for k, _ in enumerate(cat)])

                # Choose which counterpart, according to the weighting scheme
                weights = self.weight_catalogue(cat[mask])
                weights /= np.sum(weights)
                j = np.random.choice(len(weights), p=weights)
                inj_sources.append(j)

                # Simulate new source position, and remove from catalogue

                cat_obj = cat[mask][j].copy()
                cat_obj["ra_rad"], cat_obj["dec_rad"] = fixed_source.simulate_position()

                # print(fixed_source.eval_spatial_pdf(cat_obj["ra_rad"], cat_obj["dec_rad"]))
                # input("?")

                cat = np.delete(cat, j)
                inj_cat.append(cat_obj)

            inj_cat = np.array(inj_cat, dtype=cat.dtype)

            # print(np.sum(inj_cat["Flux1000"]))
            # input("?")

            cat = np.append(cat, inj_cat)

        return cat




class UniformPriorHypothesis(Hypothesis):
    name = "uniform_prior"

    @staticmethod
    def weight_catalogue(cat_data):
        return np.ones(len(cat_data))