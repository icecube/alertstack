import pickle
import os
import numpy as np
from tqdm import tqdm
from alertstack.stats import GammaDistribution
from datetime import datetime
import logging
import random
from tqdm.contrib.concurrent import process_map

import time


class Analyse:

    def __init__(self, cat, hypos, fixed_sources, cache_dir, min_E = 0, clean_cache=False):
        self.base_cat = cat
        self.fixed_sources = fixed_sources
        self.cache_dir = cache_dir
        self.hypos = dict()
        self.min_E = min_E
        for hypo in hypos:
            self.hypos[hypo.name] = hypo(fixed_sources, min_E)
        if len(hypos) == 1:
            # The fixed catalogue may change depending on the hypothesis
            self.fixed_sources = self.hypos[hypo.name].fixed_catalogue

        self._injection_hypo = None

        if clean_cache:
            self.clean_cache()

        self.all_res = dict()
        self.ts_fits = dict()
        self.sensitivity_thresholds = dict()

        self.pid = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")

    def save_path(self):
        self.pid = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
        return os.path.join(self.cache_dir, "{0}.pkl".format(self.pid))

    def set_injection_hypo(self, injection_hypo, min_E=0):
        self._injection_hypo = injection_hypo(self.fixed_sources, min_E)

    def run_trial(self, fraction=0.0, random_seed=None):
        '''
        Run individual trial. This function scrambles the catalog sources, injects the 
        correlations and calculates the TS of the trial. Returns a dictionary with the 
        results.
        '''

        # t0 = time.time()

        np.random.seed(random_seed)
        
        if fraction > 1.0:
            raise Exception("Fraction of correlated alerts cannot exceed 1.0!")

        cat = self.base_cat.scramble()
        # print(bkg_pdf_per_source)
        
        # t1 = time.time()
        # print(f"Scramble performed in: {t1-t0} s")

        if self._injection_hypo is not None:
            cat = self._injection_hypo.inject_signal(cat=cat, fraction=fraction)

        #t2 = time.time()
        #print(f"Injection performed in: {t2-t1} s")

        res = dict()

        for name, hypo in self.hypos.items():
            res[name] = [hypo.calculate_llh(
                cat, gp_threshold=self.base_cat.gp_threshold
            )]

        #t3 = time.time()
        #print(f"Llh calculated in: {t3-t2} s")

        return res

    def run_trial_wrapper(self, p):
        return self.run_trial(*p)

    def iterate_run(self, injection_hypo=None, n_trials=100, fraction=1.0, n_steps=10,
                    max_workers=min(32, os.cpu_count() + 4), chunksize=1,
                    **kwargs):
        '''
        Run the analysis. It creates the list of trials based on the input parameters. 
        It calls the run_trial function and parses the fraction of astrophysical neutrinos to inject.  
        ------------------------
        Parameters:
        :param injection_hypo: Injection Hypothesis object
        :param n_trials: Number of trials to run for each injection strength. 10x this number
        will be run as background trials
        :param fraction: Maximum fraction of astrophysical neutrinos to be injected
        :param n_steps: Number of different injection steps to test, between 0 and fraction.
        :param max_workers: tqdm max_workers parameter, setting number of cpus to be used
        :param chunksize: tqdm chunksize parameter, Size of chunks sent to worker processes
        :param kwargs: Keyword args
        '''

        self.set_injection_hypo(injection_hypo)

        # Create list of fractions to loop over. Includes ten times as many background trials.

        fs = [0.0 for _ in range(n_trials * 10)]
        for step in np.linspace(0.0, fraction, n_steps + 1)[1:]:
            fs += [step for _ in range(n_trials)]

        # Create input list

        inputs = [(x, int(random.random() * 10 ** 8)) for x in fs]

        # Run multiprocessing if circularised neutrino alerts, regular loop otherwise
 
        if 'Healpix' not in type(self.fixed_sources).__name__:
            results = process_map(self.run_trial_wrapper, inputs, max_workers=max_workers, chunksize=chunksize)
        else:
            results = []
            for i in tqdm(range(len(inputs))):
                results.append(self.run_trial(inputs[i][0],inputs[i][1])) 

        all_res = dict()

        # Combine results into nested dictionaries

        for fraction in sorted(list(set(fs))):
            mask = np.array(fs) == fraction
            cut_results = np.array(results)[mask]
            res_dict = cut_results[0]
            for entry in cut_results[1:]:
                for key, val in entry.items():
                    res_dict[key] += val

            all_res[fraction] = res_dict

        self.all_res = all_res

        self.dump_results()

    @staticmethod
    def combine_res_dicts(dict_a, dict_b):
        for hypo, hypo_res in dict_a.items():
            if hypo in dict_b.keys():
                for key, val in dict_a[hypo].items():
                    if key in dict_b[hypo].keys():
                        dict_b[hypo][key] += val
                    else:
                        dict_b[hypo][key] = val
            else:
                dict_b[hypo] = hypo_res

        return dict_b

    def dump_results(self):

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

        savepath = self.save_path()
        if os.path.isfile(savepath):
            cache_results = self.load_cache()
            self.all_res = self.combine_res_dicts(cache_results, self.all_res)

        logging.info(f"Saving to: {savepath}")

        with open(savepath, "wb") as f:
            pickle.dump(self.all_res, f)

    def load_cache(self):
        savepath = self.save_path()
        with open(savepath, "rb") as f:
            cache_results = pickle.load(f)
        return cache_results

    def find_cache_files(self):
        return [os.path.join(self.cache_dir, x) for x in os.listdir(self.cache_dir) if ".pkl" in x]

    def load_results(self, filename=None):

        self.all_res = dict()

        if filename is None:
            list_of_files = self.find_cache_files()
            latest_file = max(list_of_files, key=os.path.getctime)
        else:
            latest_file = os.path.join(self.cache_dir, filename)

        with open(latest_file, "rb") as f:
            cache_dict = pickle.load(f)
            self.all_res = self.combine_res_dicts(self.all_res, cache_dict)

        #self.clean_cache()
        self.dump_results()
        self.fit_results()
        return self.all_res

    def clean_cache(self):
        for file in self.find_cache_files():
            os.remove(file)

    def fit_results(self):
        for key, val in self.all_res[0.0].items():
            self.sensitivity_thresholds[key] = np.median(val)
            self.ts_fits[key] = GammaDistribution(val)

    def discovery_threshold(self, hypo, sigma=5.):
        return self.ts_fits[hypo].calculate_discovery_potential(sigma)

    def find_overfluctuations(self, key, threshold, **kwargs):
        pass

    def find_sensitivity(self):
        return self.find_overfluctuations("sensitivity", 0.9)

    def find_discovery_potential(self, sigma=5.):
        return self.find_overfluctuations("discovery", 0.5, sigma=sigma)







