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
from scipy.spatial.distance import pdist,squareform
from scipy.spatial import KDTree
from tqdm import tqdm
import cv2
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
        biomass (pandas.DataFrame): Pandas Dataframe containing the biomass vs time data from biomass.csv
        ntypes (pandas.DataFrame): Pandas Dataframe containing the cell number vs time data from ntypes.csv
        avg_con (pandas.DataFrame): Pandas Dataframe containing the average nutrient concentrations vs time data from avg_concentration.csv
        positions (pandas.DataFrame): Pandas Dataframe containing the single cell biomass over time of all cell ids present at the timepoint

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
            self.sucRatio = int(self.directory.split('_')[3])
        elif id:
            self.id = id
            self.get_datafed_data()
            self.sucRatio = int(self.directory.split('_')[3])
        elif directory and id:
            print('Define either a local directory or DataFed Collection ID, not both')
        else:
            print('Missing local directory or DataFed Collection ID')

        self.dims += [self.dims.pop(0)] # move Z dimension to last: z,x,y to x,y,z
        self.numsteps = len(self.timepoints)
        self.Timesteps = self.positions.Timestep.unique()
        self.calc_biomass()
        self.get_mothers()
        
    def get_local_data(self):
        """
        Collect NUFEB simulation data from a local directory.
        """
        try:
            h5 = h5py.File(os.path.join(self.directory,'trajectory.h5'),mode='r')
            self.timepoints = [key for key in h5['concentration']['co2'].keys()]
            self.timepoints.sort(key=int)
            self.dims = list(h5['concentration']['co2']['0'].shape)
            self.nutrients = list(h5['concentration'].keys())
            self.collect_positions(h5)
            self.get_nutrient_grid(h5)
            h5.close()
        except:
            print('Missing HDF5 file')

        self.biomass = pd.read_csv(os.path.join(
            self.directory,'Results','biomass.csv'),
            usecols=[0,1,2],delimiter='\t')
        self.ntypes = pd.read_csv(os.path.join(
            self.directory,'Results','ntypes.csv'),
            usecols=[0,1,2],delimiter='\t')
        self.avg_con = pd.read_csv(os.path.join(
            self.directory,'Results','avg_concentration.csv'),
            usecols=[0,2,3,4],
            delimiter='\t',
            names=['Time','O2','Sucrose','CO2'],
            skiprows=1)
        f = open(os.path.join(self.directory,'metadata.json'),'r')
        self.metadata  = json.load(f)
        f.close()
        self.convert_units_avg_con()
        self.convert_units_biomass()

    def convert_units_avg_con(self):
        """Convert the object attribute avg_con, which contains the average nutrient concentration, units to hours and mM.
        """
        self.avg_con.index = self.avg_con.Time/60/60*self.timestep
        self.avg_con.index.name='Hours'
        self.avg_con.drop('Time',inplace=True,axis=1)
        SucroseMW = 342.3
        O2MW = 32
        CO2MW = 44.01
        self.avg_con.O2 = self.avg_con.O2/O2MW*1e3
        self.avg_con.Sucrose = self.avg_con.Sucrose/SucroseMW*1e3
        self.avg_con.loc[:,'CO2'] = self.avg_con.loc[:,'CO2']/CO2MW*1e3

    def convert_units_biomass(self):
        """Convert the object attribute biomass units to hours and femtograms.
        """
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
        try:
            h5 = h5py.File(os.path.join(self.directory,'trajectory.h5'))
            self.timepoints = [key for key in h5['concentration']['co2'].keys()]
            self.timepoints.sort(key=int)
            self.dims = list(h5['concentration']['co2']['0'].shape)
            self.nutrients = list(h5['concentration'].keys())
            self.collect_positions(h5)
            self.get_nutrient_grid(h5)
            h5.close()
        except:
            print('Missing HDF5 file')

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
 

    def calc_biomass(self):
        df = self.positions
        df['biomass'] = 0
        df.loc[df.type==1,'biomass'] = 4/3*np.pi*df['radius']**3*self.metadata['cyano']['Density']*1e18
        df.loc[df.type==2,'biomass'] = 4/3*np.pi*df['radius']**3*self.metadata['ecw']['Density']*1e18
        df['time'] = df['Timestep'].transform(lambda j: j/360)

    def collect_positions(self,h5):
        """
        Extract the x, y, z position of each cell during the simulation.

        Args:
            timepoint (int):
                The simulation timestep to get the position data from.
        
        Returns:
            pandas.DataFrame:
                Dataframe containing Timestep, ID, type, radius, x, y, z columns
        """
        dfs = list()
        for t in self.timepoints:
            dfs.append(
                pd.concat([pd.Series(np.ones(h5['x'][str(t)].len())*int(t),dtype=int,name='Timestep'),
                pd.Series(h5['id'][str(t)],name='ID'),
                pd.Series(h5['type'][str(t)],name='type'),
                pd.Series(h5['radius'][str(t)],name='radius'),
                pd.Series(h5['x'][str(t)],name='x'),
                pd.Series(h5['y'][str(t)],name='y'),
                pd.Series(h5['z'][str(t)],name='z')],axis=1)
                )
        temp = pd.concat(dfs,ignore_index=True)
        idx = temp[temp.type==0].index
        self.positions = temp.drop(idx).reset_index(drop=True)


    def get_neighbor_distance(self,id,timepoint):
        """
        Get the nearest neighbor cell distances

        Args:
            id (int):
                The ID of the reference cell
            timepoint (int):
                The timepoint to check the neighbor distances from
        Returns:
            pandas.DataFrame:
                Dataframe containing ID, type, Distance
        """
        # TODO Speed up or parallelize this computation
        df = self.positions[self.positions.Timestep==timepoint]
        temp = (df.loc[df.ID ==id,['x','y','z']].squeeze() - df.loc[df.ID !=id,['x','y','z']])**2
        dist = pd.Series(np.sqrt(temp.x + temp.y + temp.z),name='Distance')
        return pd.concat([df.loc[df.ID !=id,['ID','type']],dist],axis=1).reset_index(drop=True)

    def get_neighbors(self,timestep):
        """
        Get the nearest neighbor cell distances

        Args:
            timestep (int): 
                The timepoint to check the neighbor distances from

        Returns:
            pd.DataFrame:
                Pandas dataframe containing pairwise neighbor distances
        """
        df = self.positions
        df2 = df[df.Timestep == timestep].set_index(['ID'])
        df2.sort_index(inplace=True)
        # distances =pdist(df2[['x','y','z']])
        pairwise = pd.DataFrame(
            squareform(pdist(df2[['x','y','z']])),
            columns = df2.index,
            index = df2.index
        )
        pairwise[pairwise ==0] = np.nan
        return pairwise

    def get_mothers__old(self):
        """
        Assign mother cells based on initial cells in the simulation.

        Returns:
            pandas.DataFrame:
                Dataframe containing ID, type, position, radius, and mother_cell

        """
        df = self.positions
        df['mother_cell'] = -1
        for ID in df.loc[df.Timestep==0,'ID'].unique():
            idx = df[df['ID'] ==ID].index
            df.loc[idx,'mother_cell'] = ID

        for time in tqdm(sorted(df[df.Timestep!=0].Timestep.unique()),desc='Assigning ancestry'):
            for type_ in df.type.unique():
                ancestors = df[(df.type==type_) & (df.Timestep==time) & (df.mother_cell != -1)]
                arr1 = ancestors[['x','y','z']].to_numpy()
                tree1 = KDTree(arr1)
                motherless = df[(df.type==type_) & (df.Timestep==time) & (df.mother_cell == -1)]
                if not motherless.empty:
                    d, i = tree1.query(motherless[['x','y','z']].to_numpy(), k=1)
                    idx1 =motherless.index
                    a = ancestors.iloc[i,:].mother_cell.values
                    df.loc[idx1,'mother_cell']=a
        self.colonies = df
    def get_mothers(self):
        """
        Assign mother cells based on initial cells in the simulation.

        Returns:
            pandas.DataFrame:
                Dataframe containing Timestep, ID, type, position, radius, biomass, total biomass, and mother_cell

        """
        df = self.positions.copy()
        df['mother_cell'] = -1
        df.loc[df.Timestep==0,'mother_cell'] = df.loc[df.Timestep==0,'ID']
        ancestry_df = df.loc[df.Timestep==0,['ID','mother_cell']]
        type_=1
        for time in tqdm(sorted(df[df.Timestep!=0].Timestep.unique()),desc='Assigning ancestry'):
            for type_ in df.type.unique():
                temp = df.loc[(df.type==type_) & (df.Timestep==time),['ID','x','y','z']]
                ancestors = temp.join(ancestry_df.set_index(['ID']),on='ID',how='inner', lsuffix='_left', rsuffix='_right')
                arr = ancestors[['x','y','z']].to_numpy()
                tree= KDTree(arr)
                motherless = pd.merge(temp,ancestors,on='ID',how='left', indicator=True).query('_merge == "left_only"').drop('_merge', 1).drop('x_y',1).iloc[:,:4]

                if not motherless.empty:
                    d, i = tree.query(motherless[['x_x','y_x','z_x']].to_numpy(), k=1)
                    motherless.loc[:,'mother_cell'] = ancestors.iloc[i,4].to_numpy()
                    ancestry_df = pd.concat([ancestry_df,motherless.loc[:,['ID','mother_cell']]],ignore_index=True)
        df = df.join(ancestry_df.set_index(['ID']),on='ID',how='right', lsuffix='_left', rsuffix='').drop('mother_cell_left',1)

        df['total_biomass'] = df.groupby(['mother_cell','Timestep']).cumsum()['biomass']
        self.colonies = df

    def count_colony_area(self,timestep):
        """
        Count the 2d area in pixel dimensions of each colony at a given timestep.

        Args:
            timestep (int):
                Timestep to count
        """
        if not hasattr(self,'colonies'):
            self.get_mothers()
            df = self.colonies
        else:
            df = self.colonies
        tp = df[df.Timestep == timestep]
        img_size = 2000
        bk = 255 * np.ones(shape=[img_size, img_size, 3], dtype=np.uint8)
        circles = [cv2.circle(bk,center = (round(x/self.metadata['Dimensions'][0]*img_size),
                    round(y/self.metadata['Dimensions'][1]*img_size)),radius = round(radius/self.metadata['Dimensions'][1]*img_size),
                    color = (cell,0,0),thickness = -1) for x,y, radius,cell in zip(tp.x,tp.y,tp.radius,tp.mother_cell)]
        cols, counts = np.unique(bk[:,:,0],return_counts=1)
        for colony,area in zip(cols[:-1],counts[:-1]):
            idx = df[(df.mother_cell==int(colony)) & (df.Timestep==timestep)].index
            self.colonies.loc[idx,'Colony Area'] = area

    def get_colony_areas(self):
        """Count colony areas for all timesteps
        """
        if not hasattr(self,'colonies'):
            self.get_mothers()
            df = self.colonies
        else:
            df = self.colonies
        for time in tqdm(df.Timestep.unique(),desc='Counting colony areas'):
            self.count_colony_area(time)

    def get_nutrient_grid(self,h5):
        keys = list(h5['concentration'].keys())
        timepoints = [k for k in h5['concentration'][keys[0]].keys()]
        timepoints.sort(key=int)
        stacks = list()
        for key in keys:
            dfs = list()
            for time in timepoints:
                dfs.append(h5['concentration'][key][time])
            stacks.append(np.stack(dfs))
        grid = np.stack(stacks,axis=1)
        self.grid = grid
        return
    def get_local_con(self,timestep,cellID):
        """
        Get the local nutrient concentration of a cell

        Args:

            timestep (int):
                The timestep at which to check the concentration
            cellID (int):
                The cell identification number
        
        Returns:
            Nutrient Concentration (float):
                The concentration of the specified nutrient within the cell's grid
        """
        cell_locs = self.positions
        grid = [np.linspace(0,self.metadata['Dimensions'][x],self.dims[x]) for x in range(3)]
        grid_loc = [get_grid_idx(grid[i],cell_locs[cell_locs.ID ==cellID][d].values[0]) for i,d in enumerate(['x','y','z'])]

        return self.grid[timestep,:,grid_loc[2],grid_loc[0],grid_loc[1]]

    def get_fitness(self,timestep,cellID):
        """
        Get the fitness of an individual cell based on the relative Monod growth rate at a given timestep

        Args:
            timestep (int):
                The timestep at which to check the concentration
            cellID (int):
                The cell identification number
        Returns:
            float:
                The Monod growth rate (1/s)
        """
            # TODO Speed up or parallelize this computation
        df = self.positions
        cell_type = df[(df.Timestep ==timestep) & (df.ID ==cellID)].type.values[0]
        if df[(df.Timestep ==timestep) & (df.ID ==cellID)].empty:
            print('Timestep or cell ID not found')
            return
        concentrations = self.get_local_con(list(df.Timestep.unique()).index(timestep),cellID)
        if cell_type == 1:
            metadata = self.metadata['cyano']
            light = concentrations[self.nutrients.index('sub')]
            co2 = concentrations[self.nutrients.index('co2')]
            fitness = metadata['GrowthRate'] * (light / (metadata['K_s']['sub'] + light)) * (co2 / (metadata['K_s']['co2'] + co2))
            return fitness
        elif cell_type == 2:
            metadata = self.metadata['ecw']
            suc = concentrations[self.nutrients.index('suc')]
            o2 = concentrations[self.nutrients.index('o2')]
            maintenance = metadata['GrowthParams']['Maintenance'] * (o2 / (metadata['K_s']['o2'] + o2))
            decay = metadata['GrowthParams']['Decay']
            fitness = metadata['GrowthRate'] * (suc / (metadata['K_s']['suc'] + suc)) * (o2 / (metadata['K_s']['o2'] + o2))
            return fitness - maintenance - decay
    def collect_fitness(self):
        df = self.positions
        fitness = pd.DataFrame(columns=['Time','ID','Fitness'])
        for time in tqdm(self.Timesteps):
            for cell in df[(df.Timestep==time)].ID:
                fitness = fitness.append(pd.DataFrame([[time,cell,self.get_fitness(time,cell)]],columns=['Time','ID','Fitness']),ignore_index=True)
        self.fitness=fitness

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
    from datafed.CommandLib import API
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
