import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import os

from alertstack.analyse import Analyse
from alertstack.scramble_catalogues.agn_radio_catalogue import (
    AstrogeoAGNCatalogue,
    AverageFluxWeightHypothesis,
)
from alertstack.stats import GammaDistribution, TSHandler
from examples.radio_agn_neutrino_alert import agn_analysis
from scipy import stats

cwd = os.path.dirname(os.path.abspath(__file__))
analysis_cache_dir = os.path.join(cwd, "cache/")
figures_folder = os.path.join(cwd,"figures/")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Calculate TS distributions'
    )
    parser.add_argument(
        '--n_trials',
        type=int,
        default=500,
        help = 'Number of trials'
    )
    parser.add_argument(
        '--fraction',
        type=float,
        default=0.0,
        help = 'Maximum fraction of neutrinos to be correlated'
    )
    parser.add_argument(
        '--n_steps', type=int, default=0, help ='Number of steps'
    ) 
    parser.add_argument(
        '--input',
        type=str,
        default=(
            "/data/user/gsommani/alertstack-icecube/examples/"
            "radio_agn_neutrino_alert/cache/"
            "february_update_2026_03_02-09_56_36.pkl"
        ),
        help = 'Results to use. If None, recalculates the results.'
    )
    args = parser.parse_args()
    
    '''
    n_trials: Number of trials to run for each injection strength. 10x this number
    will be run as background trials
    fraction: Maximum fraction of astrophysical neutrinos to be injected
    n_steps: Number of different injection steps to test, between 0 and fraction.
    input: Results to use. If None, recalculates the results.
    '''

    logging.getLogger().setLevel("INFO")

    # Run trials and save results
    if args.input == 'None':
        inputfile = agn_analysis.iterate_run(
            n_trials=args.n_trials,
            injection_hypo=AverageFluxWeightHypothesis,
            fraction=args.fraction,
            n_steps=args.n_steps,
            chunksize=10,
            
        )
    else:
        inputfile = args.input

    # Get the TS with the real data
    base_cat = AstrogeoAGNCatalogue()
    
    cat = base_cat.scramble() # remove this line to use the real direction of the blazars
    #cat = base_cat.data
    
    
    for name, hypo in agn_analysis.hypos.items():
        # save a list with information of correlations
        ts = hypo.calculate_llh(cat, savedata=analysis_cache_dir) 

    # Load the results from the trials
    all_res = agn_analysis.load_results(filename=inputfile)

    # Get the TS distribution of the background
    ts_handler = TSHandler(all_res, agn_analysis)
    gd = ts_handler.find_thresholds_gamma()
    key = list(ts_handler.sens_threshold.keys())[0]
    val_bkg = all_res[0][key]
    ts_handler.plot_ts(val_bkg, key, gd=gd, bins=30, ts=ts) # Plot TS distribution of bkg + TS_data

    plt.title(f"RFC2025c + IceCat-2 -> {len(
        val_bkg
    ):.1e} Scrambles")
    plt.savefig(
        figures_folder + "RFC2025c_scrambles_result",
        bbox_inches="tight",
        dpi=200
    )
    plt.savefig(
        figures_folder + "RFC2025c_scrambles_result.pdf",
        bbox_inches="tight",
        dpi=200
    )
    
    # Calculate and print p-value
    pv_int = gd.dist.sf(ts)
    pv_cnt = sum(np.array(val_bkg)[val_bkg>ts])/sum(val_bkg)
    
    print('\n###### Global p-value ################################\n')
    print(f'Integrating gamma distribution: p-value = {pv_int:.4e} ({stats.norm.ppf(1-pv_int):.2f} sigmas)') 
    print(f'Counting bins: p-value = {pv_cnt:.4e} ({stats.norm.ppf(1-pv_cnt):.2f} sigmas)')
    print('\n######################################################')
