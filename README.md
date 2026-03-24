# alertstack

## Description

Link to the wiki page: https://wiki.icecube.wisc.edu/index.php/Stacking_analysis_with_the_new_realtime_alerts_contours

This stacking analysis aims to calculate the overall correlation of neutrino alerts with blazars and tidal disruption events (TDEs). The list of neutrino alerts is fixed, using the real directions in each trial. The position of the sources to be correlated (e.g. blazars or TDEs) is randomly scrambled instead. The angular uncertainties of the neutrino alerts are the probability maps from [IceCat-2](https://wiki.icecube.wisc.edu/index.php/IceCat-2).

The main code is called [alertstack](https://github.com/icecube/alertstack) and tests for correlations between neutrino alerts and catalogs. 

## Dataset

- Neutrino catalog: [IceCat-2](https://wiki.icecube.wisc.edu/index.php/IceCat-2) (`/data/ana/realtime/alert_catalog_v3/`)
- Blazar catalog: [4LAC-DR3](https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR3/)
- AGN catalog: [Atrogeo (RFC)](http://astrogeo.org/sol/rfc/rfc_2025c/)
- Accretion flares catalog: [S. van Velzen et al. (2024)](https://academic.oup.com/mnras/article/529/3/2559/7616944?login=true)
- Infrared flares catalog: [Flaires](https://www.aanda.org/articles/aa/full_html/2025/03/aa51340-24/aa51340-24.html)

## Repository

The repository is structured as follows:

**Catalog data**
- `alertstack/data`: contains the various catalogs:
    - Several versions of the Fermi catalog:
        - 4LAC-DR1 (<code>table-4LAC.fits</code>);
        - 4LAC-DR2 (<code>table-4LAC-DR2-h.fits</code>);
        - 4LAC-DR3 (<code>table-4LAC-DR3-h.fits</code>).
    - Several versions of the Astrogeo catalog:
        - RFC 2019d (<code>rfc_2019d_cat.txt</code>);
        - RFC 2020c (<code>rfc_2020c_cat.txt</code>);
        - RFC 2022a (<code>rfc_2022a_cat.txt</code>);
        - RFC 2025b (<code>rfc_2025b_cat.txt</code>);
        - RFC 2025c, the one used in the anlysis (<code>rfc_2025c_cat.txt</code>).
    - The 63 accretion flares from [S. van Velzen et al. (2024)](https://academic.oup.com/mnras/article/529/3/2559/7616944?login=true) (<code>sjoert_catalog.txt</code> and <code>ZTF_neoWISE_flares_acflares.dat</code> for the coordinates)
    - The [Flaires catalog](https://www.aanda.org/articles/aa/full_html/2025/03/aa51340-24/aa51340-24.html) (<code>flaires.dat</code>)
    - A file with information of the light curves of blazars at the neutrino arrival times (`weights_LC.pkl`)
    - A reduced sample of public alerts used by Plavin et al. (<code>plavin_alerts.txt</code> and <code>catalog_of_alerts.txt</code>).

**General scripts**
- `alertstack/analyse.py`: defines class Analyse, which produces the background maps and runs the trials. Other useful functions include saving and loading results. 
- `alertstack/stats.py`: calculates discovery potential and sensitivity and plots all of this.
- `alertstack/plotting_tools.py`: tools to plot catalogs and data-derived distributions.

**Datasets**
- `alertstack/fixed_catalogues/icecube_neutrino_alerts.py`: loads neutrino catalog (multi-order healpix maps with probability).
- `alertstack/scrambled_catalogues/blazar_catalogue.py`: loads Fermi catalog, applies cut on latitude and energy flux.
- `alertstack/scrambled_catalogues/agn_radio_catalogue.py`: loads Astrogeo catalog, applies the cut S_8GHz > 0.15 mJy.
- `alertstack/scrambled_catalogues/sjoert_catalogue.py`: loads the catalog with the 63 accretion flares from [S. van Velzen et al. (2024)](https://academic.oup.com/mnras/article/529/3/2559/7616944?login=true) and applies same weights as in the paper.
- `alertstack/scrambled_catalogues/flaires_catalogue.py`: loads the flaires catalog, selects only the sources with an estimation of the bolometric fluence, and excludes the ones coincident with jetted synchrotron-emitting sources (normal AGN activity).

Two scripts are available to investigate the catalogs:
- `alertstack/scramble_catalogues/plot_catalogs.py`: performs all-sky scatter plots of the catalogs and shows what the data-derived background distributions look like.
- `alertstack/scramble_catalogues/ks_tests.py`: for the Astrogeo, the 63 accretion flares, and the Flaires catalog (the catalogs that make use of the data-derived background distribution), it performs a Kolmogorov-Smirnov (KS) test against uniform and data-derived distribution.

**Specific to the analyses**
- `examples/fermi_blazar_neutrino_alert`: contains the analysis scripts for the blazar-neutrino alerts correlation analysis using the average 12-year energy flux as weight.
    - `run_blazar_analysis.py`: calculates the test-statistic distribution (explained in **Run the code** section).
    - `run_blazar_analysis.sh`: wrapper of `run_blazar_analysis.py` to perform everything with one simple bash command (explained in **Run the code** section).
    - `plot_sensitivity_fermi_blazar.py`: from the test-statistic distribution, extrapolates sensitivity, discovery potential, and makes plots about them.
    - `fermi_isotropic.ipynb`: Jupyter Notebook to study the isotropy of the Fermi catalog.
   
- `examples/fermi_LC_blazar_neutrino_alert`: contains the main analysis script for the blazar-neutrino alerts correlation analysis using the available light curves as weights.
    - `run_blazar_analysis.py`: calculates the test-statistic distribution (explained in **Run the code** section).
    - `run_blazar_analysis.sh`: wrapper of `run_blazar_analysis.py` to perform everything with one simple bash command (explained in **Run the code** section).
    - `plot_sensitivity_fermi_LC_blazar.py`: from the test-statistic distribution, extrapolates sensitivity, discovery potential, and makes plots about them.
    - `create_weights_LC.py`: creates a pickle file with a dictionary containing the values of each blazar light curve at the arrival times of all neutrinos. The file is already stored in `alertstack/data`.

- `examples/radio_agn_neutrino_alert`: contains the main analysis script for the AGN-neutrino alerts correlation analysis.
    - `run_radio_analysis.py`: calculates the test-statistic distribution (explained in **Run the code** section).
    - `run_radio_analysis.sh`: wrapper of `run_radio_analysis.py` to perform everything with one simple bash command (explained in **Run the code** section).
    - `plot_sensitivity_radio_agn.py`: from the test-statistic distribution, extrapolates sensitivity, discovery potential, and makes plots about them.

- `examples/sjoert_accretion_flares_neutrino_alert`: contains the main analysis script for the 63 accretion flares-neutrino alerts correlation analysis.
    - `run_sjoert_analysis.py`: calculates the test-statistic distribution (explained in **Run the code** section).
    - `run_sjoert_analysis.sh`: wrapper of `run_sjoert_analysis.py` to perform everything with one simple bash command (explained in **Run the code** section).
    - `plot_sensitivity_sjoert.py`: from the test-statistic distribution, extrapolates sensitivity, discovery potential, and makes plots about them.
    - `dag_creator_sjoert.py`: script that creates the dag file to submit the jobs to recreate the full test statistic distribution. Since a few coincidences are enough for a very high test statistic, the test-statistic background distribution does not follow asymptotic distributions and needs to be fully explored with scrambles.

- `examples/flaires_neutrino_alert`: contains the main analysis script for the Flaires-neutrino alerts correlation analysis.
    - `run_flaires_analysis.py`: calculates the test-statistic distribution (explained in **Run the code** section).
    - `run_flaires_analysis.sh`: wrapper of `run_sjoert_analysis.py` to perform everything with one simple bash command (explained in **Run the code** section).
    - `plot_sensitivity_flaires.py`: from the test-statistic distribution, extrapolates sensitivity, discovery potential, and makes plots about them.
    - `dag_creator_flaires.py`: script that creates the dag file to submit the jobs to recreate the full test statistic distribution. Since a few coincidences are enough for a very high test statistic, the test-statistic background distribution does not follow asymptotic distributions and needs to be fully explored with scrambles.

To plot all test-statistic distributions, discover potentials, and sensitivities with a single command, the script `examples/plot_all.py` is available.

Moreover, the script `examples/compare_sensitivities.py` compares sensitivities with [the older version of this analysis](https://github.com/icecube/wg-nu-sources/tree/main/2022_Fermi_Blazars_Alerts_Stacking).


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

*alertstack* runs both with published neutrino information or internal probability skymaps. Icecube Collaboration members can run:

```
export NU_SKYMAP_DIR=/path/to/healpix/files 
```

to import and use Healpix files (on the cobalt machines, `/data/ana/realtime/alert_catalog_v3/fits/`).

## Run the code

To calculate the sensitivity and discovery potential of the analyses, do 

```
./examples/???_neutrino_alert/run_???_analysis.sh
```

This will run the analysis locally with the necessary parameters to obtain the results in the wiki page. You can also select the number of trials, the maximum fraction of neutrino alerts to correlate, and the number of steps to consider for the fraction with `--n_trials`, `--fraction`, and `--n_steps`. The code will calculate a TS distribution with `n_trials` trials (10x more for the background) for each fraction of injected astrophysical neutrinos, from 0.0 to `fraction`. For each TS distribution, the injection fraction is calculated as (step number) * `fraction`/`n_steps`.

Additionally, two further parameters are available: `--tag` and `--run`. `--tag` specifies a tag to save the results, to make them recognizable. `--run` is for testing. It specifies the last run to consider in the calculations. This can be useful to compare with older results.

With the default values (different for each analysis), the running time is of some hours. With the recommended values, if you just want to test if it works (`n_trials = 25`, `fraction = 0.25`, and `n_steps = 10`), it will take a few seconds (but the results will be statistically limited; this is just to get approximate values). With this set of parameters, the script will run 250 background trials and 25 trials for injections of 0.1, 0.2, 0.3, and 0.4 times the neutrino flux.

The script will not output anything. To investigate the distribution, sensitivity, and discovery potential, you should run (after you export the path to the neutrino maps and activate your virtual environment).

```
python examples/???_neutrino_alert/plot_sensitivity_???.py --input PATH-TO-TS-DISTRIBUTION
```

where `PATH-TO-TS-DISTRIBUTION` indicates the path to the `.pkl` file with the test-statistic distribution (usually in `examples/???_neutrino_alert/cache/NAME_FILE.pkl`, check the folder for the exact name). The default parameters in `plot_sensitivity_???.py` will reproduce the results in the wiki.

The output should also include the following lines at the end (the actual total number of neutrinos and average signalness may change from one analysis to another):

```
------- Sensitivity and discovery potential with 393 neutrino alerts (average signalness: 46.8 %) --------

Sensitivity at ? of flux, expectation of ?/393 = ?
3 Sigma discovery at ? of flux, expectation of ?/393 = ?
5 Sigma discovery at ? of flux, expectation of ?/393 = ?

----------------------------------------------------------------------------------------------------------
```

The various plots for each specific analysis will be saved in `examples/??_neutrino_alert/figures/`.


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
