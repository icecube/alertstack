import pickle
import os
import numpy as np
from tqdm import tqdm


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

        self.pid = None



    def save_path(self):
        return os.path.join(self.cache_dir, "{0}.pkl".format(self.pid))

    def run_trial(self, injection_hypo=None, fraction=0.0):

        if fraction > 1.0:
            raise Exception("Fraction of correlated alerts cannot exceed 1.0!")

        cat = self.base_cat.scramble()

        if injection_hypo is not None:
            cat = injection_hypo.inject_signal(cat=cat, fraction=fraction)

        res = dict()

        for name, hypo in self.hypos.items():
            res[name] = [hypo.calculate_llh(cat)]

        return res

    def run_trials(self, injection_hypo_class=None, n_trials=100, fraction=0.0):
        res = None
        if injection_hypo_class is not None:
            injection_hypo = injection_hypo_class(self.fixed_sources)
        else:
            injection_hypo = None
        for _ in tqdm(range(n_trials)):
            if res is None:
                res = self.run_trial(injection_hypo, fraction=fraction)
            else:
                trial_res = self.run_trial(injection_hypo, fraction=fraction)

                for key, val in trial_res.items():
                    res[key] += val

        return res

    def iterate_run(self, injection_hypo=None, n_trials=100, fraction=1.0, nsteps=10):

        steps = np.linspace(0.0, fraction, nsteps+1)[1:]

        all_res = dict()
        all_res[0.0] = self.run_trials(injection_hypo, n_trials*10, fraction=0.0)

        if np.logical_and(injection_hypo is not None, fraction > 0.0):
            for step in steps:
                all_res[step] = self.run_trials(injection_hypo, n_trials, fraction=step)

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
        savepath = self.save_path()
        if os.path.isfile(savepath):
            cache_results = self.load_cache()
            self.all_res = self.combine_res_dicts(cache_results, self.all_res)

        print("Saving to:", savepath)

        print(self.all_res.keys())

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

        for file in self.find_cache_files():
            with open(file, "rb") as f:
                cache_dict = pickle.load(f)
                self.all_res = self.combine_res_dicts(self.all_res, cache_dict)

        self.clean_cache()
        self.dump_results()
        return self.all_res

    def clean_cache(self):
        for file in self.find_cache_files():
            os.remove(file)

    def sensitivity_threshold(self):
        return






