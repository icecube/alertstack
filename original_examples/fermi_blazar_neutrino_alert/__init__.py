import os
import numpy as np
import matplotlib.pyplot as plt
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

def plot_ts(val, sens, disc_3, disc_5, gd, filename="ts_bkg.png"):
    # plot TS + gamma distribution + sensitivity and disc potential
    data = np.array(val)
    plt.hist(data, density=True, bins=50)
    x_range = np.linspace(min(data), max(data), 100)
    ylim = plt.gca().get_ylim() 
    plt.plot(x_range, gd.dist.pdf(x_range))
    plt.xlabel('TS')
    plt.axvline(sens, color = 'tab:orange', ls='--', label='Sensitivity')
    plt.axvline(disc_3, color = 'tab:orange', ls='-.', label='3sigma disc')
    plt.axvline(disc_5, color = 'tab:orange', ls=':', label='5sigma disc')
    plt.legend()
    plt.savefig(os.path.join(blazar_cache, filename),dpi=500)

def plot_ts_data(val, gd, ts_data):
    # plot TS + gamma distribution + TS of data
    data = np.array(val)
    plt.hist(data, density=True, bins=50)
    x_range = np.linspace(min(data), max(data), 100)
    ylim = plt.gca().get_ylim() 
    plt.plot(x_range, gd.dist.pdf(x_range))
    plt.yscale('log')
    plt.xlabel('TS')
    plt.axvline(ts_data, color = 'tab:orange', label='TS with 4LAC-DR2')
    plt.legend()
    plt.savefig(os.path.join(blazar_cache,'ts_data.png'),dpi=500)

