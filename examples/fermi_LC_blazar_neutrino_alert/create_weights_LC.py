from pathlib import Path
import pickle as pkl
import numpy as np
from tqdm import tqdm
import os
import json

from alertstack.scramble_catalogues.blazar_catalogue import Fermi4FGLBlazarCatalogue
from alertstack.fixed_catalogues.icecube_neutrino_alerts import HealpixNeutrinoAlertCatalogue
from alertstack import alertstack_data_dir

def get_lc(blazars):
    '''Get monthly light curves of blazars from Fermi 4LAC-DR2.

    Parameters
    ----------
    blazars: `pandas.DataFrame`
        The catalog of blazars for which the light curves are required.
    
    Returns
    mask_lc: `list[list]`
        list where each element is the mask for the correspondent blazar.
        0 if the value in values_lc is flux, 1 if it's an upper limit.
    values_lc: `list[list]`
        list where each element is the light curve for the correspondent blazar.
        average integrated flux of photons or upper limits for each time bin
    ts_lc: `list[list]`
        list where each element is another list, one for each blazar, having
        as alements arrays of 2-tuples, first element is center
        of bin time interval (in MET) and second is the likelihood TS
    
   
    If the light curve is not available, it will return None. 
    '''
    
    mypath='/data/user/gsommani/fermi_lightcurves/lightcurves'
    onlyfiles = [f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath, f))]
    justnames = []
    
    # Get 4FGL names of available light curves
    for i in onlyfiles:
        name = i[:4] + ' ' + i[7:14]
        if i[14]=='_':
            name = name + '+' + i[17:21]
        else:
            name = name + i[14:19]
        justnames.append(name)

    mask_lc = []
    values_lc = []
    ts_lc = []
    
    # Loop over blazars in 4LAC-DR2
    for b_name in blazars['Source_Name']:
        if b_name in justnames: # if the LC is available
            ind = justnames.index(b_name)
            with open(os.path.join(mypath,onlyfiles[ind])) as f:
                data = json.load(f)

            mask = []
            points = []
            
            # Take values of flux or upper limit for each time bin
            for i in data['ts']:
                if i[0] in [j[0] for j in data['flux']]:
                    # Case 1: There is a flux measurement
                    ind = [j[0] for j in data['flux']].index(i[0])
                    mask.append(0)
                    points.append(data['flux'][ind][1])
                elif i[0] in [j[0] for j in data['flux_upper_limits']]:
                    # Case 2: There is an upper limit
                    ind = [j[0] for j in data['flux_upper_limits']].index(i[0])
                    mask.append(1)
                    points.append(data['flux_upper_limits'][ind][1])
                else:
                    # Case 3: There is nothing, treat as upper limit (value not important)
                    mask.append(0)
                    points.append(0.0)

            mask_lc.append(mask) 
            values_lc.append(points) 
            ts_lc.append(data['ts'])
        else:
            mask_lc.append(None) 
            values_lc.append(None)
            ts_lc.append(None)

    return mask_lc, values_lc, ts_lc


def not_flux(c,t,b):
    '''Function to calculate weight if the neutrino arrived in a month where there is no data
    or it's just an upper limit. Checks the month before and after, if both have flux data the weight 
    is an average of both values. If only one of them is a data point, the assumed flux for the month 
    of interest is the same as that one.
    
    Parameters: 
    c: `float`
        center of the bin of neutrino arrival time
    t: `float`
        times of data points in the light curve
    b: `???`
        Element from the catalog of blazars, with lightcurves.
    '''
    
    before = c - 30
    after = c + 30

    fl_before = None
    fl_after = None

    if before in t: # points before and after of interest
        whether_flux_before = b['mask_lc'][np.where(t == before)[0][0]]
        if whether_flux_before == 0:
            fl_before = b['values_lc'][np.where(t == before)[0][0]]

    if after in t:
        whether_flux_after = b['mask_lc'][np.where(t == after)[0][0]]
        if whether_flux_after == 0:
            fl_after = b['values_lc'][np.where(t == after)[0][0]]

    if fl_before is not None and fl_after is not None:
        fl = (fl_before + fl_after)*0.5
    elif fl_before is not None:
        fl = fl_before
    elif fl_after is not None:
        fl = fl_after
    else: 
        fl = 0.0 # b['Energy_Flux100']

    return fl 


def flux_at_nu_new(
    b,
    nu_at,
):
    '''Get the value of the energy flux in the monthly time bin
    in which the neutrino arrived. Take as parameters a blazar (b)
    and the neutrino arrival time (nu_at) and returns the flux.

    Parameters
    ----------
    
    b: `???`
        Element from the catalog of blazars, with lightcurves.
    nu_at: `float`
        Neutrino arrival time.
    '''
    
    if b['values_lc'] != None:
        # Transform time information to MJD
        MJDREF = 51910 + 7.428703703703703e-4
        t = MJDREF + np.asarray([i[0] for i in b['ts_lc']])/86400
        bins = np.arange(t[0]-15, t[-1]+45, 30)
        ind = np.digitize(nu_at, bins) - 1 
        c = t[0] + ind*30
        
        if c in t: # not a gap:
            whether_flux = b['mask_lc'][np.where(t == c)[0][0]]
            if whether_flux == 0: # it's flux
                fl = b['values_lc'][np.where(t == c)[0][0]]
            else: # it's upper limit
                fl = not_flux(c,t,b)
                # value cannot be higher than upper limit or negative
                if (fl > b['values_lc'][np.where(t == c)[0][0]]) and (
                    b['values_lc'][np.where(t == c)[0][0]] > 0
                ):
                    fl = b['values_lc'][np.where(t == c)[0][0]]
        else: # in a gap
            fl = not_flux(c,t,b)
    else:
        fl = b['Energy_Flux100']
        
    return fl

# Load catalogs
blazar_cat = Fermi4FGLBlazarCatalogue()
nu_cat = HealpixNeutrinoAlertCatalogue()

# Get monthly light curves for each blazar
mask_lc, values_lc, ts_lc = get_lc(blazar_cat.data)
new_dt = np.dtype(
    blazar_cat.data.to_records(
        index=False
    ).dtype.descr + [(
        'values_lc', list
    )] + [('mask_lc', list)] + [('ts_lc', list)]
)
b = np.zeros(blazar_cat.data.to_records(index=False).shape, dtype=new_dt)
for i in blazar_cat.data.to_records(index=False).dtype.descr:
    b[i[0]] = blazar_cat.data[i[0]]

b['mask_lc'] = mask_lc
b['values_lc'] = values_lc
b['ts_lc'] = ts_lc

# Calculate weights
weights = {}
for btmp in tqdm(b):
    val = {}
    for nu in nu_cat:
        val[nu.time_mjd] = flux_at_nu_new(btmp, nu.time_mjd)
    weights[btmp['Source_Name']] = val
    
# Store values
a = Path(alertstack_data_dir) / 'weights_LC.pkl' 
with a.open('wb') as f:
    pkl.dump(weights,f)