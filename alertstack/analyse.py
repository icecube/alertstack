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
    """Class that treats the single analysis as an object.

    Parameters
    ----------
    cat: `pandas.DataFrame`
        Catalogue of astrophysical sources
    hypos: `list(alertstack.Hypothesis)`
        List of hypotheses to be tested. Initially, the list was meant
        to test more hypotheses at the same time, but now this is
        in practice never done. Maybe this should be updated?
    fixed_sources: `alertstack.FixedCatalogue`
        Catalogue of neutrino alerts
    cache_dir: `str`
        Directory where the results will be cached.
    min_E: `float`
        Cut all neutrino events below this energy (in TeV)
        (useful to investigate the minimal sensitive energy)
    clean_cache: `bool`
        If true, clean the cache directory before saving new results.
    """

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
        """Get the path where the results will be saved.
        """
        self.pid = datetime.now().strftime("%Y_%m_%d-%H_%M_%S")
        return os.path.join(self.cache_dir, "{0}.pkl".format(self.pid))

    def set_injection_hypo(self, injection_hypo, min_E=0, max_run=200000):
        """Set the hypothesis that will determine how the injections will work.
    
        Parameters
        ----------
        injection_hypo: `alertstack.Hypothesis`
            hypothesis that will determine how the injections will work
        min_E: `float`
            Cut all neutrino events below this energy (in TeV)
            (useful to investigate the minimal sensitive energy)
        max_run: `int`
            Remove all neutrino alerts after this run   
        """
        self._injection_hypo = injection_hypo(self.fixed_sources, min_E, max_run)

    def run_trial(self, fraction=0.0, random_seed=None):
        '''Run individual trial. This function scrambles the catalog sources,
        injects the correlations and calculates the TS of the trial.
        Returns a dictionary with the results.

        Parameters
        ----------
        fraction: `float`
            Fraction of astrophysical neutrino flux to inject.
        random_seed: `float | None`
            Seed used to sample randomly from distributions
        '''
        np.random.seed(random_seed)
        
        if fraction > 1.0:
            raise Exception("Fraction of correlated alerts cannot exceed 1.0!")

        cat = self.base_cat.scramble()

        if self._injection_hypo is not None:
            cat = self._injection_hypo.inject_signal(cat=cat, fraction=fraction)

        res = dict()

        for name, hypo in self.hypos.items():
            res[name] = [hypo.calculate_llh(
                cat, gp_threshold=self.base_cat.gp_threshold
            )]

        return res

    def run_trial_wrapper(self, p):
        '''Wrapper for the function run_trial. Useful to parallelize.

        Parameters
        ----------
        p: `???`
            All the parameters for the function run_trial
        '''
        return self.run_trial(*p)

    def iterate_run(
        self,
        injection_hypo=None,
        n_trials=100,
        fraction=1.0,
        n_steps=10,
        max_workers=min(32, os.cpu_count() + 4),
        chunksize=1,
        additional_tag="",
        progression_bar=True,
        max_run=200000
    ):
        '''Run the analysis. It creates the list of trials based on the input parameters. 
        It calls the run_trial function and parses the fraction of astrophysical neutrinos to inject.  
        
        Parameters:
        -----------
        injection_hypo: `alertstack.Hypothesis`
            Injection Hypothesis object
        n_trials: `int`
            Number of trials to run for each injection strength.
            10x this number will be run as background trials
        fraction: `float`
            Maximum fraction of astrophysical neutrinos to be injected
        n_steps: `int`
            Number of different injection steps to test between 0 and fraction.
        max_workers: `int`
            tqdm max_workers parameter, setting number of cpus to be used
        chunksize: `int`
            tqdm chunksize parameter, Size of chunks sent to worker processes
        additional_tag: `str`
            additional tag to add to the filename for results
        progression_bar: `bool`
            Show in real-time the progress of the iterations on the terminal
        max_run: `int`
            Remove all neutrino alerts after this run   
        '''

        self.set_injection_hypo(injection_hypo, max_run=max_run)

        # Create list of fractions to loop over. Includes ten times as many background trials.
        fs = [0.0 for _ in range(n_trials * 10)]
        for step in np.linspace(0.0, fraction, n_steps + 1)[1:]:
            fs += [step for _ in range(n_trials)]

        # Create input list
        inputs = [(x, int(random.random() * 10 ** 8)) for x in fs]

        # Run multiprocessing if circularised neutrino alerts, regular loop otherwise
        if 'Healpix' not in type(self.fixed_sources).__name__:
            if not progression_bar: 
                t0 = time.time()
            results = process_map(
                self.run_trial_wrapper,
                inputs,
                max_workers=max_workers,
                chunksize=chunksize,
                disable=not progression_bar,
            )
            if not progression_bar:
                tot_time_s = time.time() - t0
                tot_time_min = int(tot_time_s / 60)
                tot_time_hrs = int(tot_time_min / 60)
                tot_time_s = tot_time_s % 60
                print(
                    f"Iterations concluded! It took: {tot_time_hrs} h"
                    f" {tot_time_min} min {tot_time_s:.2f} s"
                )
        else:
            results = []
            if progression_bar:
                for i in tqdm(range(len(inputs))):
                    results.append(self.run_trial(inputs[i][0],inputs[i][1])) 
            else:
                t0 = time.time()
                for i in range(len(inputs)):
                    results.append(self.run_trial(inputs[i][0],inputs[i][1]))
                tot_time_s = time.time() - t0
                tot_time_min = int(tot_time_s / 60)
                tot_time_hrs = int(tot_time_min / 60)
                tot_time_s = tot_time_s % 60
                print(
                    f"Iterations concluded! It took: {tot_time_hrs} h"
                    f" {tot_time_min} min {tot_time_s:.2f} s"
                )

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

        self.dump_results(additional_tag=additional_tag)

    @staticmethod
    def combine_res_dicts(dict_a, dict_b):
        '''Combine dictionaries for different results in one single dict
        
        Parameters:
        -----------
        dict_a: `dict`
            first dictionary
        dict_b: `dict`
            second dictionary
        '''
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

    def dump_results(self, additional_tag=""):
        '''Save the results inside a file
        
        Parameters:
        -----------
        additional_tag: `str`
            Additional string to tag the specific file
        '''

        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)


        if additional_tag != "":
            additional_tag += "_"
        savepath = (
            f"{self.save_path().split(self.save_path().split("/")[-1])[0]}"
            f"{additional_tag}{self.save_path().split("/")[-1]}"
        )
        if os.path.isfile(savepath):
            cache_results = self.load_cache()
            self.all_res = self.combine_res_dicts(cache_results, self.all_res)

        logging.info(f"Saving to: {savepath}")

        with open(savepath, "wb") as f:
            pickle.dump(self.all_res, f)

    def load_cache(self):
        '''If the save_path is a file and not a directory,
        load this file.
        '''
        savepath = self.save_path()
        with open(savepath, "rb") as f:
            cache_results = pickle.load(f)
        return cache_results

    def find_cache_files(self):
        '''Find all result files in the cache directory.
        '''
        return [os.path.join(
            self.cache_dir, x
        ) for x in os.listdir(self.cache_dir) if ".pkl" in x]

    def load_results(self, filename=None, dump_results=False):
        '''Load result files.

        parameters
        ----------
        filename: `str|None`
            If none, take the latest file in the cache directory.
            If string, take the file pointed to by the path.
        dump_results: `bool`
            Dump or not the results after loading them
        '''

        self.all_res = dict()

        if filename is None:
            list_of_files = self.find_cache_files()
            latest_file = max(list_of_files, key=os.path.getctime)
        else:
            latest_file = os.path.join(self.cache_dir, filename)

        with open(latest_file, "rb") as f:
            cache_dict = pickle.load(f)
            self.all_res = self.combine_res_dicts(self.all_res, cache_dict)

        if dump_results:
            self.dump_results()
        return self.all_res

    def clean_cache(self):
        '''Clean the cache directory.
        '''
        for file in self.find_cache_files():
            os.remove(file)







