import argparse
import glob
import numpy as np
import os
import pickle

from tqdm import tqdm

analysis = "sjoert"
hypothesis_key = "strength_flux_weight"

cwd = os.path.dirname(os.path.realpath(__file__))
cache_folder = os.path.join(cwd,"cache/")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Unify the results from several jobs'
    )
    parser.add_argument(
        '--n_trials',
        type=int,
        default=100000,
        help = 'Number of trials per job'
    )
    parser.add_argument(
        '--fraction',
        type=float,
        default=0.25,
        help = 'Maximum fraction of neutrinos to be correlated'
    )
    parser.add_argument(
        '--n_steps', type=int, default=15, help ='Number of steps'
    )
    parser.add_argument(
        '--tag',
        type=str,
        default="",
        help = 'Alternative tag to identify the results'
    )
    args = parser.parse_args()

    """
    n_trials: Number of trials to run for each injection strength,
    per each job. 10x this number will be run as background trials
    fraction: Maximum fraction of astrophysical neutrinos to be injected
    n_steps: Number of different injection steps to test,
    between 0 and fraction.
    tag: Alternative tag to identify the results.
    """

    n_trials = args.n_trials
    fraction = args.fraction
    n_steps = args.n_steps
    
    if args.tag == "":
        tag = f"{analysis}_n{n_trials}_f{fraction}_s{n_steps}"
    else:
        tag = args.tag
    
    results_paths = glob.glob(os.path.join(cache_folder,f"{tag}_N*"))
    print(f"\nFound {len(results_paths)} results with the tag {tag}.\n")

    total_results = dict()
    
    for i, path in tqdm(enumerate(results_paths)):
        with open(path, 'rb') as f:
            data = pickle.load(f)
            if i == 0:
                keys = data.keys()
                res_list = [[] for j in range(len(keys))]
            for k, key in enumerate(keys):
                res_list[k] = res_list[k] + list(data[key][hypothesis_key])

    for i, key in enumerate(keys):
        total_results[key] = {
            hypothesis_key : np.array(res_list[i])
        }

    if args.tag == "":
        newtag = (
            f"{analysis}_n{n_trials*len(results_paths)}"
            f"_f{fraction}_s{n_steps}"
        )
    else:
        newtag = f"{tag}_times_{len(results_paths)}"
    filename = f"{newtag}.pkl"
    savepath = os.path.join(cache_folder, filename)
    
    with open(savepath, "wb") as f:
        pickle.dump(total_results, f)

    print(f"\nSaved the combined results in {savepath}\n")
        