import os
import argparse
import sys
from glob import glob
import h5py
from pathlib import Path
import seaborn as sns
import pandas as pd
from datafed.CommandLib import API
df_api = API()
df_api.setContext('p/eng107')

class NUFEB_data:
    """
    NUFEB simulation data class to collect results for analysis
    """
    def __init__(self,directory,local=True,id=None):
        self.directory = directory
        self.local = local
        self.id = id
        self.sucRatio = int(self.directory.split('_')[-2])
        if self.local:
            self.get_local_data()
        elif self.local is not True and self.id is not None:
            self.datafed = self.get_datafed_data(self.id)
        else:
            print('Something went wrong')
        self.timepoints = [key for key in self.h5['concentration']['co2'].keys()]
        self.timepoints.sort(key=int)
        self.dims = self.h5['concentration']['co2']['0'].shape
        self.numsteps = len(self.timepoints)
    def get_local_data(self):
        """
        Collect NUFEB simulation data from a local directory
        """
        self.h5 = h5py.File(os.path.join(self.directory,'trajectory.h5'))
        self.biomass = pd.read_csv(os.path.join(
            self.directory,'Results','biomass.csv'),
                                   usecols=[0,1,2],delimiter='\t')
        self.ntypes = pd.read_csv(os.path.join(
            self.directory,'Results','ntypes.csv'))
        self.avg_con = pd.read_csv(os.path.join(
            self.directory,'Results','avg_concentration.csv'),
            usecols=[0,2,3,4],
            delimiter='\t',
            names=['Time','O2','Sucrose','CO2'],
            skiprows=1)
    def get_datafed_data(record_id):
        """
        Collect NUFEB simulation data from a DataFed collection

        Args:
            record_id (str):
                The DataFed record identifier of the data
        """
        dv_resp = df_api.dataView(record_id)
        get_resp = df_api.dataGet(record_id,
                          '.', # directory where data should be downloaded
                          orig_fname=False, # do not name file by its original name
                          wait=True, # Wait until Globus transfer completes
                         )
        return get_resp
    #print(dv_resp)
    def radius_key(timestep):
        """
        Generate the appropriate key for a radius at a given timestep

        Args:
            timestep (str):
                The simulation timestep to get keys from

        Example:

        >>> ts='1000'
        >>> radius_key(ts)
        """
        return(f"radius{timestep}")