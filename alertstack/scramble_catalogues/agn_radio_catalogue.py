from astropy.io import fits
import copy
import healpy as hp
import numpy as np
import os
import logging
from alertstack import AnisotropicExtragalacticCatalogue, Hypothesis, is_outside_GP, alertstack_data_dir
from numpy.lib.recfunctions import rename_fields
from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS
import pandas as pd


class AstrogeoAGNCatalogue(AnisotropicExtragalacticCatalogue):
    '''Loads Astrogeo RFC catalog and selects AGNs with S > 0.15 mJy.
    '''

    @staticmethod
    def parse_data(name_cat='rfc_2025c_cat.txt'):
        """Load the catalogue.

        Parameters
        ----------        
        name_cat: `str | None`
            Possibility to specify the name of the catalog to use.
        """

        logger = logging.Logger("default_logger")
        logger.setLevel("DEBUG")

        d = []

        with open(os.path.join(alertstack_data_dir, name_cat), 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    if name_cat == 'rfc_2020c_cat.txt' or name_cat == 'rfc_2022a_cat.txt' or name_cat == 'rfc_2022b_cat.txt':
                        d.append(
                            {
                                'Category': line.split()[0], 
                                'IVS name': line.split()[1], 
                                'Source_Name': line.split()[2], 
                                'ra': [float(i) for i in line.split()[3:6]], 
                                'dec': [float(i) for i in line.split()[6:9]],  
                                'N of obs': int(line.split()[12]), 
                                'S band map': float(line.split()[13]), 
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
                    elif name_cat == 'rfc_2025c_cat.txt':
                        d.append(
                            {
                                'Comm name': line.split()[2], 
                                'Source_Name': line.split()[1], 
                                'ra': [float(i) for i in line.split()[3:6]], 
                                'dec': [float(i) for i in line.split()[6:9]],  
                                'N of obs': int(line.split()[13]), 
                                'S band map': float(line.split()[15]), 
                                'S band unresolved': line.split()[17],
                                'C band map': line.split()[18],
                                'C band unresolved': line.split()[20],            
                                'X band map': float(line.split()[21]), 
                                'X band unresolved': line.split()[23], 
                                'U band map': line.split()[24],
                                'U band unresolved': line.split()[26], 
                                'K band map': line.split()[27], 
                                'K band unresolved': line.split()[29],
                                # 'Type': line.split()[24], 
                                # 'Catalog': line.split()[25] 
                            }
                        )
            df = pd.DataFrame(d)
        

        logging.info("Selecting AGNs with S > 0.15 mJy")

        agns = df.loc[df['X band map']>=0.15]
        agns.reset_index(drop=True,inplace=True)

        agns_ras = list(agns['ra'])
        agns_decs = list(agns['dec'])
        agns_ras = ['{0:.0f}h{1:.0f}m{2}s'.format(i[0],i[1],i[2]) for i in agns_ras]
        agns_decs = ['{0:.0f}d{1:.0f}m{2}s'.format(i[0],i[1],i[2]) for i in agns_decs]
        agns_coord = [agns_ras[i] + " " + decs for i,decs in enumerate(agns_decs)]
        agns_coord = [SkyCoord(i, frame=ICRS) for i in agns_coord]
        new_decs = np.array([c.dec.deg for c in agns_coord])
        new_ras = np.array([c.ra.deg for c in agns_coord])
        
        maps = {
            "ra": "ra_deg",
            "dec": "dec_deg",
        }
        
        agns = agns.rename(columns=maps)
        agns["ra_deg"] = new_ras
        agns["dec_deg"] = new_decs
        agns.insert(len(agns.keys()), 'ra_rad', new_ras * np.pi / 180.)
        agns.insert(len(agns.keys()), 'dec_rad', new_decs * np.pi / 180.) # called like this but it's in deg actually (same for 4LAC)
        agns.insert(
            len(agns.keys()), 'bkg_pdf', np.empty(len(agns))
        )

        logging.info("Found {0} sources in total".format(len(agns)))
        
        # agns = agns.to_records(index = False)

        return agns

    @staticmethod
    def set_gp_threshold():
        """Set a cut in galactic latitude for the catalogue
        (exclude the sources with a smaller latitude in absolute value).
        """
        return 0.

    @staticmethod
    def set_min_declination():
        """Set a cut in declination for the catalogue
        (exclude the sources with a smaller declination).
        """
        return -90.

    @staticmethod
    def set_nside():
        """Set the nside for the final resolution of the healpix map
        describing the distribution of sources.
        """
        return 128

    def set_bkg_distribution(self):
        """Contains the logic necessary to generate an appropriate
        background distribution for the specific catalogue.
        """
        return  self.apply_cuts_on_bkg_distribution(
            self.generate_bkg_distribution_allsky(
                self.data, nside=16, hd_nside=self.nside, sigma_smoothing=8.,
            ),
            self.min_declination,
            self.gp_threshold,
        )


class AverageFluxWeightHypothesis(Hypothesis):
    """Class for the hypothesis of constant emission proportional
    to the X band flux.
    """
    name = "average_radio_flux_weight"
    unit = "mJy"

    @staticmethod
    def weight_catalogue(cat_data):
        """Weight the astrophysical sources according to the hypothesis.

        Parameters
        ----------
        cat_data: `pandas.DataFrame`
            catalogue to weight
        """
        return cat_data['X band map']