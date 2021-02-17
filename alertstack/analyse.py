import pickle
import os
import numpy as np
from tqdm import tqdm
from alertstack.stats import GammaDistribution
from datetime import datetime
import logging
import random
from tqdm.contrib.concurrent import process_map


class Analyse:

    def __init__(self, cat, hypos, fixed_sources, cache_dir, clean_cache=False):
        self.base_cat = cat
        self.fixed_sources = fixed_sources
        self.cache_dir = cache_dir
        self.hypos = dict()
        for hypo in hypos:
            self.hypos[hypo.name] = hypo(fixed_sources)

        if clean_cache:
            self.clean_cache()

        self.all_res = dict()
        self.ts_fits = dict()
        self.sensitivity_thresholds = dict()

        self.pid = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")

    def save_path(self):
        self.pid = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
        return os.path.join(self.cache_dir, "{0}.pkl".format(self.pid))

    def run_trial(self, injection_hypo=None, fraction=0.0, random_seed=None):

        np.random.seed(random_seed)

        if fraction > 1.0:
            raise Exception("Fraction of correlated alerts cannot exceed 1.0!")

        cat = self.base_cat.scramble()

        if injection_hypo is not None:
            cat = injection_hypo.inject_signal(cat=cat, fraction=fraction)

        res = dict()

        for name, hypo in self.hypos.items():
            res[name] = [hypo.calculate_llh(cat)]

        return res

    def run_trial_wrapper(self, p):
        return self.run_trial(*p)

    def iterate_run(self, injection_hypo=None, n_trials=100, fraction=1.0, n_steps=10, **kwargs):

        fs = [0.0 for _ in range(n_trials * 10)]
        for step in np.linspace(0.0, fraction, n_steps + 1)[1:]:
            fs += [step for _ in range(n_trials)]

        inputs = [(injection_hypo, x, int(random.random() * 10 ** 8)) for x in fs]
        results = process_map(self.run_trial_wrapper, inputs, **kwargs)
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

    def load_results(self):

        self.all_res = dict()

        list_of_files = self.find_cache_files()
        latest_file = max(list_of_files, key=os.path.getctime)
        #for file in self.find_cache_files():
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







