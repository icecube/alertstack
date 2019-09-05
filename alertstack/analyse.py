import healpy as hp
import random
import numpy as np
from tqdm import tqdm


class Analyse:

    def __init__(self, cat, hypos, fixed_sources):
        self.base_cat = cat
        self.fixed_sources = fixed_sources
        self.hypos = dict()
        for hypo in hypos:
            self.hypos[hypo.name] = hypo(fixed_sources)

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

        return all_res





