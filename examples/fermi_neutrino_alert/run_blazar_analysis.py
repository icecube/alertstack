from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue, AverageFluxWeightHypothesis
from alertstack.fixed_catalogues.icecube_neutrino_alerts import CircularisedNeutrinoAlertCatalogue

ana = Analyse(
    Fermi4FGLBlazarCatalogue(),
    [AverageFluxWeightHypothesis],
    CircularisedNeutrinoAlertCatalogue()
)

ana.run_trial()