import argparse
import logging

from alertstack.scramble_catalogues.flaires_catalogue import (
    FluencebolHypothesis
)
from examples.flaires_neutrino_alert import (
    flaires_analysis
)

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
        default=0.5,
        help = 'Maximum fraction of neutrinos to be correlated'
    )
    parser.add_argument(
        '--n_steps', type=int, default=10, help ='Number of steps'
    )
    parser.add_argument(
        '--tag', type=str, default=None, help ='Additional tag'
    )
    parser.add_argument(
        '--run',
        type=int,
        default=200000,
        help ='Last run to consider for the neutrinos'
    )
    parser.add_argument(
        '--progression_bar',
        type=bool,
        default=False,
        help ='Show or not the progression bar'
    )
    args = parser.parse_args()
    
    '''
    n_trials: Number of trials to run for each injection strength. 10x this number
    will be run as background trials
    fraction: Maximum fraction of astrophysical neutrinos to be injected
    n_steps: Number of different injection steps to test, between 0 and fraction.
    run: Last run to consider for the neutrinos.
    progression_bar: Show or not the progression bar
    '''

    logging.getLogger().setLevel("INFO")

    chunksize = 10

    # Run analysis and save results
    flaires_analysis.iterate_run(
        n_trials=args.n_trials,
        injection_hypo=FluencebolHypothesis,
        fraction=args.fraction,
        n_steps=args.n_steps,
        additional_tag=args.tag,
        chunksize=chunksize,
        progression_bar=args.progression_bar
    )