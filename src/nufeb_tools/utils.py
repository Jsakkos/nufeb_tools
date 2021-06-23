import os
import h5py
from pathlib import Path
import pandas as pd
import numpy as np
import json
import subprocess
import sys
import argparse
import pickle
import json
from urllib.parse import urlparse
from urllib.request import urlretrieve
import tarfile
from nufeb_tools import __version__

urls= ['https://github.com/Jsakkos/nufeb-tools/raw/main/data/runs.tar']

class get_data:
    """Collect results for analysis.

    NUFEB simulation data class to collect results for analysis

    Attributes:
        test (bool): Set `test = True` to get example data from the Github repository
        directory (str): Path to the directory containing NUFEB simulation data. If `directory = None`, get_data will look for a DataFed collection
        id (str):DataFed record ID, e.g., `"c/34558900"`
        timestep (int): Length of simulation timestep in seconds
        SucRatio (int): Relative cyanobacterial sucrose secretion level, 0-100
        timepoints (List(str)): List of timepoints in the simulation
        dims (List(str)): Size of the simulation boundaries in micrometers
        numsteps (int): Number of timepoints
        h5 (HDF5 File): HDF5 file containing NUFEB simulation data
        biomass (pandas.DataFrame): Pandas Dataframe containing the biomass vs time data from biomass.csv
        ntypes (pandas.DataFrame): Pandas Dataframe containing the cell number vs time data from ntypes.csv
        avg_con (pandas.DataFrame): Pandas Dataframe containing the average nutrient concentrations vs time data from avg_concentration.csv
        single_cell_biomass (pandas.DataFrame): Pandas Dataframe containing the single cell biomass over time of all cell ids present at the timepoint

    """
    def __init__(self,directory=None,id=None,test=None,timestep=10):
        self.timestep=timestep
        if test:
            self.directory = str((Path.home()) / '.nufeb_tools' / 'data' / 'Run_60_18_63_1_2021-06-23')
            if not os.path.isdir(self.directory):
                download_test_data()     
            self.get_local_data()
        elif directory:
            self.directory = directory
            self.get_local_data()
            self.sucRatio = int(self.directory.split('_')[-2])
        elif id:
            self.id = id
            self.get_datafed_data()
            self.sucRatio = int(self.directory.split('_')[3])
        elif directory and id:
            print('Define either a local directory or DataFed Collection ID, not both')
        else:
            print('Missing local directory or DataFed Collection ID')
        try:
            self.timepoints = [key for key in self.h5['concentration']['co2'].keys()]
            self.timepoints.sort(key=int)
            self.dims = list(self.h5['concentration']['co2']['0'].shape)
            self.dims += [self.dims.pop(0)] # move Z dimension to last: z,x,y to x,y,z
            self.numsteps = len(self.timepoints)
        except AttributeError:
            print('Missing HDF5 file')
        
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
        f = open(os.path.join(self.directory,'metadata.json'),)
        self.metadata  = json.load(f)
        f.close()
        self.convert_units_avg_con()
        self.convert_units_biomass()
    def convert_units_avg_con(self):
        self.avg_con.index = self.avg_con.Time/60/60*self.timestep
        self.avg_con.index.name='Hours'
        self.avg_con.drop('Time',inplace=True,axis=1)
        SucroseMW = 342.3
        O2MW = 32
        CO2MW = 44.01
        self.avg_con.O2 = self.avg_con.O2/O2MW*1e3
        self.avg_con.Sucrose = self.avg_con.Sucrose/SucroseMW*1e3
        self.avg_con['CO2'] = self.avg_con['CO2']/CO2MW*1e3
    def convert_units_biomass(self):
        self.biomass.index = self.biomass.step/60/60*self.timestep
        self.biomass.index.name='Hours'
        self.biomass.iloc[:,1:]=self.biomass.iloc[:,1:]*1e18
    def get_datafed_data(self,dir=None,orig_fname=True,wait=True):
        """
        Collect NUFEB simulation data from a DataFed collection

        Args:
            dir (str):
                Directory to download the DataFed collection to. Defaults to user/.nufeb/data/collection
            orig_fname (bool):
                Use original filenames
            wait (bool):
                Wait for the download to complete before moving on
        """
        from datafed.CommandLib import API
        df_api = API()
        df_api.setContext('p/eng107')
        if not dir:
            dv_resp = df_api.collectionView(self.id)
            self.directory = str((Path.home()) / '.nufeb_tools' / 'data' / dv_resp[0].coll[0].title)
        get_resp = df_api.dataGet(self.id,
                            path=self.directory,
                            orig_fname=orig_fname,
                            wait=wait,
                            )
        self.h5 = h5py.File(os.path.join(self.directory,'trajectory.h5'))
        self.metadata = json.loads(dv_resp[0].data[0].metadata)
        self.biomass = pd.read_csv(os.path.join(
            self.directory,'biomass.csv'),
                                   usecols=[0,1,2],delimiter='\t')
        self.ntypes = pd.read_csv(os.path.join(
            self.directory,'ntypes.csv'))
        self.avg_con = pd.read_csv(os.path.join(
            self.directory,'avg_concentration.csv'),
            usecols=[0,2,3,4],
            delimiter='\t',
            names=['Time','O2','Sucrose','CO2'],
            skiprows=1)
        self.convert_units_avg_con()
        self.convert_units_biomass()
    #print(dv_resp)
    def radius_key(self,timestep):
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
    def single_cell_growth(self,timepoint=0):
        """
        Extract single cell biomass over time from the HDF5 data. 

        Args:
            timepoint (int):
                The simulation timestep to calculate from. Default = 0.
        """
        self.hrs = [int(x)/360 for x in self.timepoints]
        df = pd.DataFrame(columns=['id','type','time','biomass'])
        # loop over all cell ids, c, from time = 0
        for c in self.h5['id'][str(timepoint)]:
            # loop over all timepoints t
            for t,h in zip(self.timepoints,self.hrs):
                # get list of ids from timepoint t
                ids = self.h5['id'][t]
                arr = ids.__array__(ids.dtype)
                # make sure cell id matches id of interest
                i = np.where(arr == c)[0][0]
                # get radius
                radius = self.h5[self.radius_key(t)][i]
                volume = 4/3*np.pi*radius**3 #volume in m^3
                # get celltype
                celltype=self.h5['type'][t][i]
                # color cells
                if celltype==1:
                    color = '#2ca25f'
                    mass = volume*370*1e18 # convert mass in kg to fg
                elif celltype ==2:
                    color = '#de2d26'
                    mass = volume*236*1e18
                elif celltype ==0:
                    print('Celltype is 0',i,c,celltype)
                # append data to a dataframe
                df = df.append(pd.DataFrame([[c,celltype,h,mass]],columns=['id','type','time','biomass']),ignore_index=True)
        self.single_cell_biomass = df
    def get_positions(self,timepoint=0):
        """
        Extract the x, y, z position of each cell at a given timepoint

        Args:
            timepoint (int):
                The simulation timestep to get the position data from.
        
        Returns:
            df (pandas.DataFrame):
                Dataframe containing ID, type,x, y, z columns
        """
        return pd.concat([pd.Series(self.h5['id'][str(timepoint)],name='ID'),
        pd.Series(self.h5['type'][str(timepoint)],name='type'),
        pd.Series(self.h5['x'][str(timepoint)],name='x'),
        pd.Series(self.h5['y'][str(timepoint)],name='y'),
        pd.Series(self.h5['z'][str(timepoint)],name='z')],axis=1)
    def get_neighbor_distance(self,id,timepoint):
        """
        Get the nearest neighbor cell distances

        Args:
            id (int):
                The ID of the reference cell
            timepoint (int):
                The timepoint to check the neighbor distances from
        Returns:
            df (pandas.DataFrame):
                Dataframe containing ID, type, Distance
        """
        df = self.get_positions(timepoint)
        temp = (df[df.ID ==id][['x','y','z']].squeeze() - df[df.ID !=id][['x','y','z']])**2
        dist = pd.Series(np.sqrt(temp.x + temp.y + temp.z),name='Distance')
        return pd.concat([df[df.ID !=id][['ID','type']],dist],axis=1).reset_index(drop=True)

    def get_local_con(self,nutrient,timestep,cellID):
        """
        Get the local nutrient concentration of a cell

        Args:
            nutrient (str):
                The nutrient to check the concentration of
            timestep (int):
                The timestep at which to check the concentration
            cellID (int):
                The cell identification number
        
        Returns:
            Nutrient Concentration (float):
                The concentration of the specified nutrient within the cell's grid
        """
        cell_locs = self.get_positions(timestep)
        grid = [np.linspace(0,self.metadata['Dimensions'][x],self.dims[x]) for x in range(3)]
        grid_loc = [get_grid_idx(grid[i],cell_locs[cell_locs.ID ==cellID][d].values[0]) for i,d in enumerate(['x','y','z'])]
        if self.h5['concentration'].__contains__(nutrient):
            return self.h5['concentration'][nutrient][str(timestep)][grid_loc[2],grid_loc[0],grid_loc[1]]
        else:
            nutes = list(self.h5['concentration'].__iter__())
            print('Nutrient ' + nutrient + ' not found. Try:')
            print(*nutes)
    def get_fitness(self,timestep,cellID):
        """
        Get the fitness of an individual cell based on the relative Monod growth rate at a given timestep

        Args:
            timestep (int):
                The timestep at which to check the concentration
            cellID (int):
                The cell identification number
        Returns:
            fitness (float):
                The Monod growth rate (1/s)
        """
        if self.h5['type'].__contains__(str(timestep)): 
            cell_type = self.h5['type'][str(timestep)][np.where(self.h5['id'][str(timestep)].__array__() == cellID)[0][0]]
        else:
            print('Timestep or cell ID not found')
            return
        if cell_type == 1:
            metadata = self.metadata['cyano']
            light = self.get_local_con('sub',timestep,cellID)
            co2 = self.get_local_con('co2',timestep,cellID)
            fitness = metadata['GrowthRate'] * (light / (metadata['K_s']['sub'] + light)) * (co2 / (metadata['K_s']['co2'] + co2))
            return fitness
        elif cell_type == 2:
            metadata = self.metadata['ecw']
            suc = self.get_local_con('suc',timestep,cellID)
            o2 = self.get_local_con('o2',timestep,cellID)
            maintenance = metadata['GrowthParams']['Maintenance'] * (o2 / (metadata['K_s']['o2'] + o2))
            decay = metadata['GrowthParams']['Decay']
            fitness = metadata['GrowthRate'] * (suc / (metadata['K_s']['suc'] + suc)) * (o2 / (metadata['K_s']['o2'] + o2))
            return fitness - maintenance - decay

def get_grid_idx(array,value):
    """
    Find the nutrient grid index value. Taken from https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array.

    Args:
        array (numpy.array):
            1D Array containing the grid positions
        value (float):
            Cell location to map to the grid
    Returns:
        index (int):
            Grid index
    """
    n = len(array)

    jl = 0# Initialize lower
    ju = n-1# and upper limits.
    while (ju-jl > 1):# If we are not yet done,
        jm=(ju+jl) >> 1# compute a midpoint with a bitshift
        if (value >= array[jm]):
            jl=jm# and replace either the lower limit
        else:
            ju=jm# or the upper limit, as appropriate.
        # Repeat until the test condition is satisfied.
    if (value == array[0]):# edge cases at bottom
        return 0
    elif (value == array[n-1]):# and top
        return n-1
    else:
        return jl             
def download_test_data(urls=urls):
    """
    Get an example dataset from the Github repo. Downloads to "home/.nufeb_tools/data"

    Args:
        urls (List(str))
    """
    # nufeb_tools directory
    cp_dir = Path.home().joinpath('.nufeb_tools')
    cp_dir.mkdir(exist_ok=True)
    data_dir = cp_dir.joinpath('data')
    data_dir.mkdir(exist_ok=True)
    # TODO Add progress bar
    for url in urls:
        parts = urlparse(url)
        filename = os.path.basename(parts.path)
        cached_file = os.path.join(data_dir, filename)
        if not os.path.exists(cached_file):
            local_filename, headers = urlretrieve(url, cached_file)
            tar = tarfile.open(local_filename,'r')
            tar.extractall(path=data_dir)
            tar.close()
            Path(local_filename).unlink()



def upload_datafed(file, title, collection_id,metadata_file):
    """
    Create a data collection to hold NUFEB data in DataFed

    Args:
        file (str):
            Path of file to upload
        title (str):
            Name to use for file on DataFed
        collection_id (str):
            The identifier of the collection to store the file
        metadata_file (str):
            Path of the metadata file to append to the data file
    """
    filename = file
    file_title= title
    global_coll_id = collection_id
    df_api = API()
    if metadata_file != '':
        pkl_file = metadata_file

        with open(pkl_file, 'rb') as f:
            metadata = pickle.load(f)
        rec_msg = df_api.dataCreate(title = file_title,
                                    alias = '',
                                    metadata=json.dumps(metadata),
                                    parent_id=global_coll_id,
                                        )
        rec_id = rec_msg[0].data[0].id
        #Use as pathname the path and name of the file you wish to move from CADES to DataFed
        pput_msg = df_api.dataPut(rec_id, filename, wait=False)
        #_logger.info(pput_msg)
    else:
        #_logger.debug('No metadata file found')
        sys.exit(1)

def create_datafed_collection(n_cyanos, n_ecw, SucPct,dims):
    """
    Create a data collection to hold NUFEB data in DataFed

    Args:
        n_cyanos (int): 
            Number of initial cyanobacteria
        n_ecw (int):
            Number of initial E. coli
        SucPct (int):
            Percentage of sucrose secretion activation
        dims (List(float)):
            x, y, z simulation boundaries
    """
    try:
        from datafed.CommandLib import API
        df_api = API()
        df_api.setContext('p/eng107')
        collectionName = f'NUFEB_{n_cyanos}_{n_ecw}_{SucPct}_{dims[0]}_{dims[1]}_{dims[2]}'
        parent_collection = df_api.getAuthUser().split('/')[1]
        coll_msg = df_api.collectionCreate(collectionName,
                                        parent_id=parent_collection)
        global_coll_id = coll_msg[0].coll[0].id
        #_logger.info(global_coll_id)
    except:
        global_coll_id = None
        #_logger.debug('Unable to create collection')
    return global_coll_id

def verify_datafed_connection():
    """
    Verify Datafed installation and connection
    """
    
    try:
        from datafed.CommandLib import API
    except ImportError:
       # _logger.info('datafed not found. Installing from pip.')
        subprocess.call([sys.executable, "-m", "pip", "install", 'datafed'])
        from datafed.CommandLib import API
    df_api = API()
        #print('Success! You have DataFed: ' + df_ver)
    # Verify user authentication
    if not df_api.getAuthUser():
       print('You have not authenticated into DataFed Client')


    # Check default Globus endpoint
    if not df_api.endpointDefaultGet():
        endpoint = 'cades#CADES-OR'
        df_api.endpointDefaultSet(endpoint)



    #print('Your default Globus Endpoint in DataFed is:\n' + df_api.endpointDefaultGet())
    # Test the endpoint
    path = str((Path.home()) / '.nufeb_tools' / 'datafed')
    cp_dir = (Path.home()) / '.nufeb_tools' / 'datafed'
    cp_dir.mkdir(exist_ok=True)
    dget_resp = df_api.dataGet('d/35437908',
                            path,
                            wait=True)
    #  _logger.debug(dget_resp)
    if dget_resp[0].task[0].status == 3:
        file = (Path.home()) / '.nufeb_tools' / 'datafed' /'35437908.md5sum'
        file.unlink()
    else:
        if dget_resp[0].task[0].msg == "globus connect offline":
            print('You need to activate your Globus Endpoint and/or ensure Globus Connect Personal is running.\n'
                'Please visit https://globus.org to activate your Endpoint')
            sys.exit(1)
        elif dget_resp[0].task[0].msg == "permission denied":
            print('Globus does not have write access to this directory. \n'
                'If you are using Globus Connect Personal, ensure that this notebook runs within'
                'one of the directories where Globus has write access. You may consider moving this'
                'notebook to a valid directory or add this directory to the Globus Connect Personal settings')
            sys.exit(1)
        else:
            NotImplementedError('Get in touch with us or consider looking online to find a solution to this problem:\n' + dget_resp[0].task[0].msg)
            sys.exit(1)