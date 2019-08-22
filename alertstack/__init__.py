import numpy as np
import healpy as hp
from scipy.stats import norm


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
        on the surface of a shpere between the points
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

class Catalogue:

    def __init__(self):
        self.data = self.parse_data()

    @staticmethod
    def parse_data():
        return NotImplementedError


class FixedCatalogue(Catalogue):
    pass



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
        self.nside = nside
        print("NSIDE = {0}, Max Pixel Radius = {1} deg".format(nside, np.degrees(hp.max_pixrad(nside))))
        self.cone_ids = [x for x in range(hp.nside2npix(self.nside)) if x not in gal_plane_1024]

    def scramble_positions(self):
        indexes = np.random.choice(len(self.cone_ids), size=len(self.data))
        return self.extract_ra_dec(self.nside, indexes)

class Hypothesis:
    name = None

    def __init__(self, fixed_catalogue):
        self.fixed_catalogue = fixed_catalogue
        self.source_weights = np.array([source.eval_source_weight() for source in fixed_catalogue])
        self.source_weights /= np.mean(self.source_weights)

    @staticmethod
    def weight_catalogue(cat_data):
        return NotImplementedError

    def calculate_llh(self, cat_data):
        cat_weights = self.weight_catalogue(cat_data)
        cat_weights /= np.sum(cat_weights)

        # print(max(cat_weights))

        llh = 0

        for i, source in enumerate(self.fixed_catalogue):
            spatial_pdf = source.eval_spatial_pdf(cat_data["ra_rad"], cat_data["dec_rad"])
            source_weight = self.source_weights[i]
            llh += np.log(np.sum(source_weight * spatial_pdf * cat_weights))
            print(np.log(np.sum(source_weight * spatial_pdf * cat_weights)))
        print("LLH: {0}".format(llh))
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

        n_exp = fraction * float(len(self.fixed_catalogue))
        n_inj = np.random.poisson(n_exp)

        if n_inj > 0:

            # Choose which fixed source will have a counterpart

            # Choose which counterpart, according to the weighting scheme
            weights = self.weight_catalogue(cat)

            # Inject the counterpart based on the spatial PDF





class UniformPriorHypothesis(Hypothesis):
    name = "uniform_prior"

    @staticmethod
    def weight_catalogue(cat_data):
        return np.ones(len(cat_data))