# alertstack

## Description

Link to the wiki page: https://wiki.icecube.wisc.edu/index.php/Correlation_of_blazars_and_neutrino_alerts
This stacking analysis aims to calculate the overall correlation between neutrino alerts and blazars. The list of neutrino alerts is fixed, using the real directions in each trial. The position of the sources to be correlated (e.g. blazars) is randomly scrambled instead. The angular uncertainties of the neutrino alerts are the best-guess Millipede errors.

The main code is called [alertstack](https://github.com/icecube/alertstack) and tests for correlations between neutrino alerts and catalogs. 

## Dataset

- Neutrino catalog: v2 alert catalog (`/data/ana/realtime/alert_catalog_v2/`)
- Blazar catalog: [4LAC-DR2](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR2/)
- AGN catalog: [Atrogeo (RFC)](http://astrogeo.org/sol/rfc/rfc_2022a/)

## Repository

The repository is structured as follows:

**General scripts**
- `alertstack/data`: contains the Fermi catalog (4LAC-DR, `table-4LAC-DR2-h.fits`), the Astrogeo catalog (RFC 2022a, `rfc_2022a_cat.txt`), a file with information of the light curves of blazars at the neutrino arrival times (`weights_LC.pkl`) and a reduced sample of public alerts to run with the toy model.
- `alertstack/analyse.py`: defines class Analyse, that produces the background maps and runs the trials. Other useful functions such as saving and loading results. 
- `alertstack/stats.py`: calculates discovery potential adjusting the data to a gamma distribution.

**Datasets**
- `alertstack/fixed_catalogues/icecube_neutrino_alerts.py`: loads neutrino catalog (healpix maps with likelihood or circularized errors for the toy model).
- `alertstack/scrambled_catalogues/blazar_catalogue.py`: loads Fermi catalog, applies cut on latitude and energy flux.
- `alertstack/scrambled_catalogues/agn_radio_catalogue.py`: loads Astrogeo catalog, applies the cut S_8GHz > 0.15 mJy.

**Specific to the analyses**
- `examples/fermi_blazar_neutrino_alert`: contains the analysis scripts for the blazar-neutrino alerts correlation analysis using the average 10-year energy flux as weight.
    - `run_blazar_analysis.py`: calculates sensitivity and discovery potential (explained in **Run the code** section).
    - `unblind_blazar_analysis.py`: calculates the p-value (only run after the unblinding has been approved).
    - `create_table_correlations.py`: prints a table with all the correlations found in the data (run and event number and alert name for the neutrino, J2000 name for the source and contribution to the TS) (only run after the unblinding has been approved).
    - `blazar_toy_example.ipynb`
   
- `examples/fermi_LC_blazar_neutrino_alert`: contains the main analysis script for the blazar-neutrino alerts correlation analysis using the available light curves as weights.
    - `run_blazar_analysis.py`: calculates sensitivity and discovery potential (explained in **Run the code** section).
    - `unblind_blazar_analysis.py`: calculates the p-value (only run after the unblinding has been approved).
    - `create_table_correlations.py`: prints a table with all the correlations found in the data (run and event number and alert name for the neutrino, J2000 name for the source and contribution to the TS) (only run after the unblinding has been approved).
    - `blazar_toy_example.ipynb`
    - `create_weights_LC.py`: creates a pickle file with a dictionary containing the values of each blazar light curve at the arrival times of all neutrinos. The file is already stored in `alertstack/data`.

- `examples/radio_agn_neutrino_alert`: contains the main analysis script for the AGN-neutrino alerts correlation analysis.
    - `run_radio_agn_analysis.py`: calculates sensitivity and discovery potential (explained in **Run the code** section).
    - `unblind_radio_agn_analysis.py`: calculates the p-value (only run after the unblinding has been approved).
    - `create_table_correlations.py`: prints a table with all the correlations found in the data (run and event number and alert name for the neutrino, J2000 name for the source and contribution to the TS) (only run after the unblinding has been approved).
    - `mimic_plavin_paper.py`: runs background trials and calculates the p-value using the method from Plavin et al. paper to create the error region of the neutrinos (explained in **Run the code** section)(only run after the unblinding has been approved). 




## Install locally

The requirements to run the code are the following:
```
    python >= 3.7,
    numpy >= 1.17.0,
    healpy,
    mhealpy,
    scipy,
    matplotlib,
    astropy,
    pandas,
    coveralls,
    tqdm >= 4.42.0,
    "hellolancel @ git+https://github.com/sjoertvv/HelloLancel.git@main",
```

All of the necessary packages will be installed if you run the following command to install *alertstack* locally:

```
python3.12 -m venv alertstack-venv
source alertstack-venv/bin/activate
pip install "numpy<2"
pip install -U pip setuptools wheel
pip install https://github.com/pschella/k3match/archive/51a49a83d36bd5289bcd1c03296cf20531b4c924.zip --no-build-isolation
pip install -e alertstack/ --no-build-isolation
```

## Use full skymaps

*alertstack* runs both with published neutrino information or internal likelihood skymaps. Icecube Collaboration members can run:

```
export NU_SKYMAP_DIR=/path/to/healpix/files 
```

to import and use Healpix files (on the cobalt machines, `/data/ana/analyses/NuSources/2022_Fermi_Blazars_Alerts_Stacking`).

## Run the code

To calculate the sensitivity and discovery potential of the analyses, do 

```
python examples/???_neutrino_alert/run_???_analysis.py
```

This will run the analysis locally with the necessary parameters to obtain the results in the wiki page. You can also select the number of trials, the maximum fraction of neutrino alerts to have a correlation and the amount of steps to consider for the fraction with `--n_trials`, `--fraction` and `--n_steps`. The code will calculate a TS distribution with `n_trials` trials (10x more for the background) for each fraction of injected astrophysical neutrinos, from 0.0 to `fraction`. For each TS distribution, the injection fraction is calculated as (step number) * `fraction`/`n_steps`.

With the default values (`n_trials = 500`, `fraction = 0.5` and `n_steps = 10`) the running time is ~8hr. With the recommended values if you just want to test if it works (`n_trials = 25`, `fraction = 0.4` and `n_steps = 4`) it will take ~30min (but the results will be statistically limited, this is just to get approximate values). Taking the lattest set of parameters, the script will run 250 background trials and 25 trials for injections of 0.1, 0.2, 0.3 and 0.4 times the neutrino flux.

*To run the code on NPX* you have to request 26GB of memory (this is being investigated and the code will be improved in the near future so that less memory is needed).

The output would include the following lines at the end:

```
------- Sensitivity and discovery potential with 275 neutrino alerts (average signalness: 45.1 %) --------

Sensitivity at ? of flux, expectation of ?/275 = ?
3 Sigma discovery at ? of flux, expectation of ?/275 = ?
5 Sigma discovery at ? of flux, expectation of ?/275 = ?

----------------------------------------------------------------------------------------------------------
```

A plot of the TS distribution of the background and the sensitivity and discovery potentials will be saved in `examples/??_neutrino_alert/cache/ts_bkg.png`.

In `examples/fermi_blazar_neutrino_alert/blazar_toy_example.ipynb` you can run a toy version of the code using a limited sample of circularized public errors instead of the likelihood maps.


When you are ready to unblind the analysis, you can run

```
python examples/???_neutrino_alert/unblind_???_analysis.py
```

and the output will be 

```
###### Global p-value ################################

Integrating gamma distribution: p-value = ?? (?? sigmas)
Counting bins: p-value = ?? (?? sigmas)

######################################################
```

For the unblinding script there is an extra argument that you can parse, `--no_run`. This will make the script look for the latest file in `/cache` and recover the background trials from there. If you don't use this flag, the script will run 5000 background trials by default (which can also be modified with the tag `--n_trials`). A plot that shows the TS distribution of the background and the TS value of the data will be stored in `examples/???_neutrino_alert/cache/ts_data.png`. Moreover, a file with all the correlations found will be stored in `examples/???_neutrino_alert/cache/correlations.pkl`. With this, you can then run 


```
python examples/???_neutrino_alert/create_table_correlations.py
```

and a list with all the correlations will be printed.

For the radio catalog analysis there is an extra script, `mimic_plavin_paper.py`. In that script, the TS is calculated using Plavin's method to obtain `S_spatial` (more info in the wiki). When you run the script with 

```
python examples/agn_radio_neutrino_alert/mimic_plavin_paper.py
```

the script produces background trials and calculates the TS of the data. The output is the p-value and a plot with the TS distribution of the background and the TS of the data is stored in `examples/agn_radio_neutrino_alert/cache/ts_data_plavin.png`. You can select the number of background trials with `--n_trials` (default is 5000). 
