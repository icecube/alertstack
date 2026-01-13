import scipy
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

class GammaDistribution:
    '''
    This class receives the background TS distribution and fits it to a gamma distribution.
    The function 'calculate_discovery_potential' calculates the value of TS needed to get the discovery potential
    '''

    def __init__(self, data):

        # prepare the data
        default_loc = min(data) - 1.
        data = np.array(data)
        cut = min(data)
        mask = data > cut

        self.frac_under = np.sum(~mask)/float(len(mask))
        N_left = np.sum(~mask)

        # define initial guess parameters and bounds for minimization
        p_start = [9., default_loc, 0.5]
        p_bounds = [(0, None),
                    (None, default_loc + 0.99),
                    (1e-5, 1e5)
                    ]

        # define likelihood function to minimize to obtain the parameters that 
        # better match the data
        def func(p):
            dist = scipy.stats.gamma(p[0], loc=p[1], scale=p[2])
            loglh = dist.logpdf(data).sum()
            loglh += N_left * dist.cdf(cut)
            return -loglh

        # minimize 
        self.res = scipy.optimize.minimize(func, x0=p_start, bounds=p_bounds)
        print(self.res)
        # define gamma distribution that represents the background TS distribution
        self.dist = scipy.stats.gamma(self.res["x"][0], loc=self.res["x"][1], scale=self.res["x"][2])

    def calculate_discovery_potential(self, sigma=5.):
        threshold = (norm.cdf(sigma) - self.frac_under)/(1 - self.frac_under)
        return self.dist.ppf(threshold)