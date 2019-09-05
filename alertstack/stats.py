import scipy
from scipy.stats import norm
import numpy as np
import matplotlib.pyplot as plt

class GammaDistribution:

    def __init__(self, data):

        # default_loc = min(data) - 1e-9
        default_loc = min(data) - 1.
        data = np.array(data)
        cut = min(data)
        mask = data > cut


        self.frac_under = np.sum(~mask)/float(len(mask))
        # self.frac_under = 0.

        # left = data > min(data)

        N_left = np.sum(~mask)
        # data = data[mask]
        weights = np.ones_like(data) / float(len(data))
        mask = np.ones_like(data, dtype=np.int)

        p_start = [9., default_loc, 0.5]
        p_bounds = [(0, None),
                    (None, default_loc + 0.99),
                    (1e-5, 1e5)
                    ]

        def func(p):
            dist = scipy.stats.gamma(p[0], loc=p[1], scale=p[2])
            loglh = dist.logpdf(data).sum()
            loglh += N_left * dist.cdf(cut)
            return -loglh

        self.res = scipy.optimize.minimize(func, x0=p_start, bounds=p_bounds)
        print(self.res)
        self.dist = scipy.stats.gamma(self.res["x"][0], loc=self.res["x"][1], scale=self.res["x"][2])

    def calculate_discovery_potential(self, sigma=5.):
        threshold = (norm.cdf(sigma) - self.frac_under)/(1 - self.frac_under)
        return self.dist.ppf(threshold)

# def plot_background_ts_distribution(ts_array, path, ts_type="Standard",
#                                     ts_val=None):
#
#     try:
#         os.makedirs(os.path.dirname(path))
#     except OSError:
#         pass
#
#     ts_array = np.array(ts_array)
#     ts_array = ts_array[~np.isnan(ts_array)]
#
#     if np.sum(np.isnan(ts_array)) > 0:
#         print("TS distribution has", np.sum(np.isnan(ts_array)), "nan entries.")
#
#     fig = plt.figure()
#
#     df, loc, scale, frac_over = fit_background_ts(ts_array, ts_type)
#
#     frac_under = 1 - frac_over
#
#     five_sigma = (raw_five_sigma - frac_under) / (1. - frac_under)
#
#     plt.axhline(frac_over * (1 - five_sigma), color="r", linestyle="--")
#
#     max_ts = np.max(ts_array)
#
#     disc_potential = scipy.stats.chi2.ppf(five_sigma, df, loc, scale)
#
#     x_range = np.linspace(0., max(max_ts, disc_potential), 100)
#
#     plt.plot(x_range, frac_over * scipy.stats.chi2.pdf(x_range, df, loc, scale),
#              color="blue", label=r"$\chi^{2}$ Distribution")
#
#     def integral(x):
#
#         return (frac_under * np.sign(x) + frac_over *
#                 (scipy.stats.chi2.cdf(x, df, loc, scale)))
#
#     plt.plot(x_range, 1. - integral(x_range), color="green", linestyle="--",
#              label=r"1 - $\int f(x)$ (p-value)")
#
#     plt.axvline(disc_potential, color="r", label=r"5 $\sigma$ Threshold")
#
#     if ts_val is not None:
#         print("\n")
#
#         if not isinstance(ts_val, float):
#             ts_val = float(ts_val[0])
#
#         # print
#
#         print("Quantifying TS:", "{:.2f}".format(ts_val))
#
#         if ts_val > np.median(ts_array):
#
#             val = (ts_val - frac_under) / (1. - frac_under)
#
#             cdf = frac_under + frac_over * scipy.stats.chi2.cdf(
#                 val, df, loc, scale)
#
#             sig = norm.ppf(cdf)
#
#         else:
#             cdf = 0.
#             sig = 0.
#
#         print("Pre-trial P-value is", "{:.2E}".format(1-cdf), 1-cdf)
#         print("Significance is", "{:.2f}".format(sig), "Sigma")
#         print("\n")
#
#         plt.axvline(ts_val, color="purple",
#                     label="{:.2f}".format(ts_val) + " TS/" +
#                     "{:.2f}".format(sig) + r" $\sigma$")
#
#     else:
#         plt.annotate(
#             '{:.1f}'.format(100 * frac_under) + "% of data in delta. \n" +
#             r"$\chi^{2}$ Distribution:" + "\n   * d.o.f.=" + \
#             '{:.2f}'.format(df) + ",\n  * loc=" + '{:.2f}'.format(loc) + \
#             " \n * scale=" + '{:.2f}'.format(scale),
#             xy=(0.1, 0.2), xycoords="axes fraction", fontsize=8)
#
#     yrange = min(1. / (float(len(ts_array)) * n_bins),
#                  scipy.stats.chi2.pdf(disc_potential, df, loc, scale))
#
#     plt.yscale("log")
#     plt.xlabel(r"Test Statistic ($\lambda$)")
#     plt.legend(loc="upper right")
#     plt.savefig(path)
#     plt.close()
#
#     return disc_potential