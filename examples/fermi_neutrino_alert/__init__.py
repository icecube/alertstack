import os
from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis,\
    BrightestFluxWeightHypothesis
from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue, HealpixNeutrinoAlertCatalogue

blazar_cache = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cache/")

blazar_analysis = Analyse(
    Fermi4FGLBlazarCatalogue(),
    [AverageFluxWeightHypothesis],
    HealpixNeutrinoAlertCatalogue(),
    cache_dir=blazar_cache,
    clean_cache=False
)
