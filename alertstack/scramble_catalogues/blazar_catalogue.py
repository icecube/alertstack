from astropy.io import fits
import numpy as np
import os
import logging
from alertstack import IsotropicExtragalacticCatalogue, Hypothesis, is_outside_GP, alertstack_data_dir
from numpy.lib.recfunctions import rename_fields
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS

class Fermi4FGLBlazarCatalogue(IsotropicExtragalacticCatalogue):

    @staticmethod
    def parse_data():

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        with fits.open(os.path.join(alertstack_data_dir, "table-4LAC-DR2-h.fits")) as hdul: # table_4LAC
            cat = hdul[1].data
        cat = np.sort(cat, order="Energy_Flux100")[::-1]#Energy_Flux100

        logging.info("Selecting blazars from 4FGL catalogue")

        blazar_class = ["bll", "BLL", "fsrq", "FSRQ", "bcu", "BCU"]

        logging.info("Using all sources from class {0}".format(blazar_class))
        #
        mask = np.array([x["CLASS"] in blazar_class for x in cat])
        blazars = np.array(cat[mask])

        cut_e = -11.6
        mask_e = np.array([x["Energy_Flux100"]>10**cut_e for x in blazars])
        blazars = np.array(blazars[mask_e])

        maps = [
            ("RAJ2000", "ra_rad"),
            ("DEJ2000", "dec_rad"),
        ]

        for (old_key, new_key) in maps:

            blazars = rename_fields(blazars, {old_key: new_key})

        mask_GP = [is_outside_GP(blazars['ra_rad'][i],blazars['dec_rad'][i]) for i in range(len(blazars))]
        blazars = blazars[mask_GP] 
        
        ####################
        #fluxes = blazars['Energy_Flux100'].byteswap().newbyteorder()
        #blazar_df = pd.DataFrame({'flux': fluxes})
        #blazar_df['bins'] = pd.cut(blazar_df['flux'],50)
        #blazar_df2 = blazar_df.groupby('bins').agg({'flux': sum, 'bins': 'count'}).rename(columns = {'bins': 'count', 'flux':'background'}).reset_index()
        #blazar_df2['flux2'] = (len(blazars)/sum(blazar_df2['background']))*blazar_df2['background']
        #blazar_df2['s/b'] = blazar_df2['flux2']/(blazar_df2['count'])
        #blazar_df = blazar_df.merge(blazar_df2,how='left',on='bins')
        #blazars = np.lib.recfunctions.append_fields(blazars, 's/b', blazar_df['s/b'].values)

        logging.info("Found {0} sources in total".format(len(blazars)))

        return blazars

    def scramble(self):
        ra, dec = self.scramble_positions()
        cat = np.copy(self.data)
        cat['ra_rad'] = ra
        cat["dec_rad"] = dec
        return cat

class AstrogeoBlazarCatalogue(IsotropicExtragalacticCatalogue):

    @staticmethod
    def parse_data():

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        d = []
        with open(os.path.join(alertstack_data_dir,'rfc_2019d_cat.txt'), 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    d.append(
                        {
                            'Category': line.split()[0], 
                            'IVS name': line.split()[1], 
                            'J2000 name': line.split()[2], 
                            'ra': [float(i) for i in line.split()[3:6]], 
                            'dec': [float(i) for i in line.split()[6:9]],  
                            'N of obs': int(line.split()[12]), 
                            'S band map': line.split()[13], 
                            'S band unresolved': line.split()[14],
                            'C band map': line.split()[15],
                            'C band unresolved': line.split()[16],            
                            'X band map': float(line.split()[17]), 
                            'X band unresolved': line.split()[18], 
                            'U band map': line.split()[19],
                            'U band unresolved': line.split()[20], 
                            'K band map': line.split()[21], 
                            'K band unresolved': line.split()[22],
                            'Type': line.split()[23], 
                            'Catalog': line.split()[24] 
                        }
                    )
            df = pd.DataFrame(d)
        
        #cat = np.sort(cat, order="Energy_Flux100")[::-1] # to do

        logging.info("Selecting blazars with S > 0.15 mJy")

        blazars = df.loc[df['X band map']>=0.15]
        blazars.reset_index(drop=True,inplace=True)

        blazars_ras = list(blazars['ra'])
        blazars_decs = list(blazars['dec'])
        blazars_ras = ['{0:.0f}h{1:.0f}m{2}s'.format(i[0],i[1],i[2]) for i in blazars_ras]
        blazars_decs = ['{0:.0f}d{1:.0f}m{2}s'.format(i[0],i[1],i[2]) for i in blazars_decs]
        blazars_coord = [blazars_ras[i] + " " + decs for i,decs in enumerate(blazars_decs)]
        blazars_coord = [SkyCoord(i, frame=ICRS) for i in blazars_coord]
        blazars['ra_rad'] = [c.ra.deg for c in blazars_coord] # called like this but it's in deg actually (same for 4LAC)
        blazars['dec_rad'] = [c.dec.deg for c in blazars_coord]

        mask_GP = [is_outside_GP(blazars['ra_rad'][i],blazars['dec_rad'][i]) for i in range(len(blazars))]
        blazars = blazars[mask_GP] 

        logging.info("Found {0} sources in total".format(len(blazars)))
        
        blazars = blazars.to_records(index = False)

        return blazars
        
    def scramble(self):
        ra, dec = self.scramble_positions()
        cat = np.copy(self.data)
        #cat = self.data.copy()
        cat['ra_rad'] = ra
        cat['dec_rad'] = dec
        return cat

class AverageFluxWeightHypothesis(Hypothesis):
    name = "average_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        try:
            return cat_data["Energy_Flux100"]#Energy_Flux100
        except:
            return cat_data['X band map']

class BrightestFluxWeightHypothesis(Hypothesis):
    name = "brightest_flux_weight"

    @staticmethod
    def weight_catalogue(cat_data):
        weights = cat_data["Energy_Flux100"]
        weights[100:] = 0.
        return weights

