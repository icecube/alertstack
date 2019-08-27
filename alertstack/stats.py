import scipy
import numpy as np

class Chi2(object):

    """ A class similar to the ones from scipy.stats
       allowing to fit left-truncated chi^2 distributions.
    """

    def __init__(self, data):
        """ Fit the given ensemble of measurements with a chi^2 function.
        `data` is a list of test statistics values.
        `cut` defines where the distribution is truncated.
        """

        # data = np.mean(data)

        # three parameters will be fitted: dof, location, scale
        p_start = [1., np.mean(data), 1.]
        p_bounds = [(0., None),  # dof > 0
                    (0., max(data)),  # location < 0 for 'truncated'
                    # effect
                    (1e-5, 1e5)]  # shape ~ free

        # define the fit function: likelihood for chi^2 distribution,
        def func(p):
            dist = scipy.stats.chi2(p[0], loc=p[1], scale=p[2])
            loglh = dist.logpdf(data).sum()
            return -loglh

        res = scipy.optimize.minimize(func, x0=p_start, bounds=p_bounds)
        print(res)

        if not res.success:
            print('Chi2 fit did not converge! Result is likely garbage.')

        # self._q_left = N_left / float(N_all)
        self._cut = 0.
        self._f = scipy.stats.chi2(res.x[0], loc=min(data), scale=res.x[1])
        self._ks = scipy.stats.kstest(data, self._f.cdf)[0]
        self.ndof = res.x[0]
        self.loc = min(data)
        self.scale = res.x[1]