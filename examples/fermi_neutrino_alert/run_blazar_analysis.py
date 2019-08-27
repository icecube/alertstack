import numpy as np
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis
from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue
from alertstack.stats import Chi2

ana = Analyse(
    Fermi4FGLBlazarCatalogue(),
    [AverageFluxWeightHypothesis],
    CircularisedNeutrinoAlertCatalogue()
)

all_res = ana.iterate_run(injection_hypo=AverageFluxWeightHypothesis, fraction=0.2, nsteps=3)

sens_threshold = dict()
disc_threshold = dict()

zero_key = 0.0

for key, val in all_res[zero_key].items():
    sens_threshold[key] = np.median(val)
    # print(Chi2(val))
    # input("?")


for step, res in all_res.items():
    print("Fraction of neutrino alerts correlated to source: {0} \n".format(step))

    bkgs = dict()

    for key, val in res.items():
        print(key, np.mean(val), np.median(val), np.std(val))
        val = np.array(val)
        print("Fraction above median background:", np.sum(val > sens_threshold[key])/float(len(val)))