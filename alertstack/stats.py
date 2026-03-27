import matplotlib.pyplot as plt
import numpy as np
import scipy

from scipy import interpolate
from scipy.optimize import bisect
from scipy.stats import norm


class TSHandler:
    '''Class to handle the TS results from the various
    analyses.

    Parameters
    ----------
    results: `dict`
        Results from an analysis
    analysis: `alertstack.Analyse`
        The analysis involved for the statistics
    '''

    PROB_3S = 1.35e-3
    PROB_5S = 2.87e-7

    def __init__(self, results, analysis, max_run=142135):
        self.results = results
        self.sens_threshold = dict()
        self.disc_3_threshold = dict()
        self.disc_5_threshold = dict()
        runs = np.array([nu.runid for nu in analysis.fixed_sources])
        tmp = np.array(
            [i.weight for i in analysis.fixed_sources]
        )[runs<=max_run]
        self.avg_signalness = np.mean(tmp)
        self.n_events = len(tmp)
        self.x1 = None
        self.x2 = None
        self.x3 = None
        self.fracs = []
        self.sens = []
        self.sig3 = []
        self.sig5 = []


    def find_thresholds_gamma(self):
        '''Find thresholds in TS for sensitivity, 3 sigma, and 5 sigma
        discovery potential assuming a gamma distribution.
        '''
        
        for key, val in self.results[0.].items():
            self.sens_threshold[key] = np.median(val)
            gd = GammaDistribution(val)
            self.disc_3_threshold[key] = gd.calculate_discovery_potential(3.)
            self.disc_5_threshold[key] = gd.calculate_discovery_potential(5.)

        return gd

    def find_thresholds_from_data(self):
        '''Find thresholds in TS for sensitivity, 3 sigma, and 5 sigma
        discovery potential using the data and not assuming any
        distribution.
        '''
        
        for key, val in self.results[0.].items():
            self.sens_threshold[key] = np.median(val)
            decreasing_vals = np.flip(np.sort(val))
            done_3s = False
            done_5s = False
            if len(val) <= int(1./self.PROB_3S):
                print("Not enough scrambles for evaluating a 3 and 5 sigma level")
                self.disc_3_threshold[key] = decreasing_vals[0]
                self.disc_5_threshold[key] = decreasing_vals[0]
                done_3s = True
                done_5s = True
            elif len(val) <= int(1./self.PROB_5S):
                print("Not enough scrambles for evaluating a 5 sigma level")
                self.disc_5_threshold[key] = decreasing_vals[0]
                done_5s = True
        
            for i, v in enumerate(decreasing_vals):
                prob = (i + 1) / len(val)
                if prob >= self.PROB_3S and not done_3s:
                    self.disc_3_threshold[key] = v
                    done_3s = True
                if prob >= self.PROB_5S and not done_5s:
                    self.disc_5_threshold[key] = v
                    done_5s = True

    
    def plot_ts(self, val, key, gd=None, bins=30, density=True, ts=None):
        '''plot TS + sensitivity and disc potential

        Parameters
        ----------
        val: `list`
            The TS data
        key: `str`
            Key that describes the investigated model
        gd: `None | alertstack.stats.GammaDistribution`
            If given, gamma distribution that fits the TS
        bins: `int`
            Number of bins for the histogram
        density: `bool`
            Shows density or absolute number of scrambles
        ts: `float | None`
            real test statistic to plot
        '''
        sens = self.sens_threshold[key]
        disc_3 = self.disc_3_threshold[key]
        disc_5 = self.disc_5_threshold[key]
        data = np.array(val)
        counts, bins_pos, _ = plt.hist(data, bins=bins, alpha=0.)
        bins_centers = (bins_pos[:-1] + bins_pos[1:])/2.
        bins_widths = (bins_pos[1:] - bins_pos[:-1])/2.
        if density:
            db = np.array(np.diff(bins_pos), float)
            counts_sum = counts.sum()
            counts =  counts / db / counts_sum
            errs = np.sqrt(counts * db * counts_sum) / (db * counts_sum)
            plt.ylabel("Density")
        else:
            errs = np.sqrt(counts)
            plt.ylabel("Scrambles")
        x_range = np.logspace(np.log10(min(data[data>0.])), np.log10(max(data)), 100)
        if gd is not None:
            x_range = np.logspace(
                np.log10(min(data[data>0.])),
                np.log10(max(list(data)+[disc_5])),
                100
            )
            plt.plot(x_range, gd.dist.pdf(x_range))
            plt.ylim(gd.dist.pdf(disc_5)/4, max(counts)*4)
        if ts is not None:
            plt.axvline(ts, color="red", linewidth=2, label="Real data")
        plt.errorbar(
            bins_centers,
            counts,
            errs,
            bins_widths,
            linestyle="",
            color="black",
        )
        plt.ylim(min(counts[counts!=0.])/8, max(counts[counts!=0.])*4)
        plt.xlabel('TS')
        plt.yscale('log')
        plt.axvline(sens, color = 'tab:orange', ls='--', label='Sensitivity')
        plt.axvline(disc_3, color = 'tab:orange', ls='-.', label='3sigma disc')
        plt.axvline(disc_5, color = 'tab:orange', ls='dotted', label='5sigma disc')
        plt.legend()


    def extract_sens_dp(self, extent=0.2):
        """Extrapolate percentages of astrophysical neutrino flux
        necessary to get a TS higher than signalness, 3-sigma,
        and 5 sigma discovery potential.

        Parameters
        ----------
        extent: `float`
            max extent to extrapolate percentages.
        """
        
        levels = [
            ("Background Median", self.sens_threshold),
            ("3 Sigma Discovery Potential", self.disc_3_threshold),
            ("5 Sigma Discovery Potential", self.disc_5_threshold)
        ]
        
        above = dict()
        
        # Calculate for each step the fraction of trials that are above
        # a certain threshold and store the information
        # in the 'above' dictionary
        for step, res in self.results.items():
            print(
                "\nFraction of neutrino alerts correlated"
                " to source: {0} \n".format(step)
            )
        
            bkgs = dict()
            temp = []
        
            for key, val in res.items():
                val = np.array(val)
        
                for name, thresh in levels:
                    print("Fraction above {0}: {1}".format(
                        name, np.sum(val > thresh[key])/float(len(val))))
                    temp.append(np.sum(val > thresh[key])/float(len(val)))        
            above[step] = temp
        
        self.fracs = list(above.keys())
        self.sens = [list(
            above.values()
        )[i][0] for i in range(len(self.fracs))]
        self.sig3 = [list(
            above.values()
        )[i][1] for i in range(len(self.fracs))]
        self.sig5 = [list(
            above.values()
        )[i][2] for i in range(len(self.fracs))]
        
        # Interpolate a curve to the data points of
        # the fraction of trials above each threshold
        f1 = interpolate.interp1d(self.fracs, self.sens, kind='cubic')
        f2 = interpolate.interp1d(self.fracs, self.sig3, kind='cubic')
        f3 = interpolate.interp1d(self.fracs, self.sig5, kind='cubic')
        
        # Calculate flux needed to achieve sensitivity
        # and discovery potential
        print(
            "\n------- Sensitivity and discovery potential with "
            "{0} neutrino alerts"
            " (average signalness: {1:.1f} %) --------\n".format(
                self.n_events, 100*self.avg_signalness
            )
        )
        
        self.x1 = bisect(lambda x: f1(x)-0.9, 0, extent, xtol=1e-6)
        self.x2 = bisect(lambda x: f2(x)-0.5, 0, extent, xtol=1e-6)
        self.x3 = bisect(lambda x: f3(x)-0.5, 0, extent, xtol=1e-6)
        
        print(
            "Sensitivity at {0:.3f} of flux,"
            " expectation of {1:.1f}/{2} = {3:.3f}".format(
                self.x1,
                (self.x1*self.avg_signalness)*self.n_events,
                self.n_events,
                (self.x1*self.avg_signalness))
        )
        print(
            "3 Sigma discovery at {0:.3f} of flux,"
            " expectation of {1:.1f}/{2} = {3:.3f}".format(
                self.x2,
                (self.x2*self.avg_signalness)*self.n_events,
                self.n_events,
                (self.x2*self.avg_signalness)
            )
        )
        print(
            "5 Sigma discovery at {0:.3f} of flux,"
            " expectation of {1:.1f}/{2} = {3:.3f}".format(
                self.x3,
                (self.x3*self.avg_signalness)*self.n_events,
                self.n_events,
                (self.x3*self.avg_signalness)
            )
        )
        print(
            "\n------------------------------------------------------"
            "----------------------------------------------------\n"
        )

    
    def plot_sens_dp(self, data_derived=True):
        """Plot results of signal injections to get percentages of
        astrophysical neutrino flux for test statistic higher
        than signalness, 3-sigma, and 5-sigma discovery potential

        Parameters
        ----------
        data_derived: `bool`
            Specifies if the background distribution is data-derived or
            assumed as isotropic
        """
        
        plt.axhline(0.5, color="black", linestyle="dotted")
        plt.axhline(0.9, color="black", linestyle="dashed")
        
        plt.axvline(
            self.x1,
            linestyle="dashdot",
            label=f"Sensitivity ({(
                self.x1*self.avg_signalness
            )*self.n_events:.1f} alerts)"
        )
        plt.axvline(
            self.x2,
            linestyle="dashdot",
            color="tab:orange",
            label=r"3-$\sigma$ discovery potential"
            f"\n({(self.x2*self.avg_signalness)*self.n_events:.1f} alerts)"
        )
        plt.axvline(
            self.x3,
            linestyle="dashdot",
            color="tab:green",
            label=r"5-$\sigma$ discovery potential"
            f"\n({(self.x3*self.avg_signalness)*self.n_events:.1f} alerts)"
        )
        data_derived_string = "Assumption of isotropic catalog"
        if data_derived:
            data_derived_string = "Data-derived background distribution"
        plt.title(
            f"{self.n_events} neutrino alerts\n"
            f"{data_derived_string}\n"
            f"average signalness {round(100*self.avg_signalness,1)}%, "
            f"~{round(self.n_events*self.avg_signalness)}"
            " astrophysical neutrinos"
        )
        plt.plot(
            self.fracs,
            self.sens,
            marker='o',
            label="Fraction above median TS"
        )
        plt.plot(
            self.fracs,
            self.sig3,
            marker='o',
            label=r"Fraction above 3-$\sigma$ level"
        )
        plt.plot(
            self.fracs,
            self.sig5,
            marker='o',
            label=r"Fraction above 5-$\sigma$ level"
        )
        plt.ylabel("Fraction of samples")
        plt.xlabel(
            "Fraction of astrophysical neutrino alerts correlated to source"
        )


class GammaDistribution:
    '''This class receives the background TS distribution
    and fits it to a gamma distribution. The function
    'calculate_discovery_potential' calculates the value
    of TS needed to get the discovery potential

    Parameters
    ----------
    data: `numpy.array | list`
        The array with the TS values
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
        self.dist = scipy.stats.gamma(
            self.res["x"][0],
            loc=self.res["x"][1],
            scale=self.res["x"][2]
        )

    def calculate_discovery_potential(self, sigma=5.):
        """Calculate the discovery given the Gamma distribution

        Parameters
        ----------
        sigma: `float`
            Number of sigmas corresponding to the discovery potential.
        """
        threshold = (norm.cdf(sigma) - self.frac_under)/(1 - self.frac_under)
        return self.dist.ppf(threshold)