import argparse
import logging

from alertstack.scramble_catalogues.blazar_catalogue import (
    MonthlyFluxWeightHypothesis
)
from examples.fermi_LC_blazar_neutrino_alert import (
    blazar_analysis
)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Calculate TS distributions')
    parser.add_argument('--n_trials', type=int, default=20000,
                            help = 'Number of trials')
    parser.add_argument('--fraction', type=float, default=0.175,
                            help = 'Maximum fraction of neutrinos to be correlated')
    parser.add_argument('--n_steps', type=int, default=10, help ='Number of steps')
    parser.add_argument('--tag', type=str, default="", help ='Additional tag')
    args = parser.parse_args()
    
    '''
    n_trials: Number of trials to run for each injection strength. 10x this number
    will be run as background trials
    fraction: Maximum fraction of astrophysical neutrinos to be injected
    n_steps: Number of different injection steps to test, between 0 and fraction.
    '''

    logging.getLogger().setLevel("INFO")

    chunksize = 10

    # Run analysis and save results
    blazar_analysis.iterate_run(
        n_trials=args.n_trials,
        injection_hypo=MonthlyFluxWeightHypothesis,
        fraction=args.fraction,
        n_steps=args.n_steps,
        additional_tag=args.tag,
        chunksize=chunksize,
        progression_bar=True
    )