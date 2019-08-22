import healpy as hp
import random
import numpy as np



class Analyse:

    def __init__(self, cat, hypos, fixed_sources):
        self.base_cat = cat
        self.fixed_sources = fixed_sources
        self.hypos = dict()
        for hypo in hypos:
            self.hypos[hypo.name] = hypo(fixed_sources)

    def run_trial(self, injection_hypo=None, fraction=0.):

        cat = self.base_cat.scramble()

        if injection_hypo is not None:
            cat = injection_hypo.inject(cat, fraction)

        res = dict()

        print("Running trial!")

        for name, hypo in self.hypos.items():
            res[name] = hypo.calculate_llh(cat)





