
from scipy.spatial import KDTree
from scipy.spatial import Voronoi
import numpy as np
import pandas as pd
from nufeb_tools import __version__
# TODO modify to account for edge effects
def fitness_metrics(obj):
    """
    Function to calculate colony-level fitness metrics. 

    Mother cell: Cell ID which seeded the colony

    Type: Cell type - cyanobacteria are type 1  and E. coli are type 2

    Voronoi area: Area of species-specific Voronoi Tesselation at the beginning of the simulation

    IPTG: Sucrose induction level

    total biomass: Biomass of each colony at the end of the simulation (fg)

    Nearest 1: Distance to nearest cyanobacteria colony

    Nearest 2: Distance to nearest E. coli colony 

    Nearest Neighbor: Distance to nearest colony

    IC1: Average distance to nearest cyanobacteria colony

    IC2: Average distance to nearest E. coli colony

    IC: Average distance to nearest colony

    Relative Neighbor Dist 1: Distance to nearest cyanobacteria colony divided by IC1

    Relative Neighbor Dist 2: Distance to nearest E. coli colony divided by IC2

    Relative Neighbor Dist: Distance to nearest colony divided by IC

    Z1: Relative neighbor distance 1 divided by sqrt(D_sucrose/mu_cyano)

    Z2: Relative neighbor distance 2 divided by sqrt(D_sucrose/mu_ecw)

    Z1_2: Relative neighbor distance 1 divided by sqrt(D_sucrose/mu_ecw)

    Z2_1: Relative neighbor distance 2 divided by sqrt(D_sucrose/mu_cyano)

    LogNearest 1: log(Nearest 1)

    LogNearest 2: log(Nearest 2)

    LogNearest: log(Neatest Neighbot)

    Inv1: Inverse sum of neighbor distance 1

    Inv2: Inverse sum of neighbor distance 2

    Log Inv1: Log squared inverse sum of neighbor distance 1

    Log Inv2: Log squared inverse sum of neighbor distance 2
    
    Colony Area: 2D area of colony at the end of the simulation

    Args:
        obj (nufeb_tools.utils.get_data): 
            Data object collected with nufeb_tools.utils.get_data

    Returns:
        pandas.DataFrame: 
            Dataframe containing colony number (mother cell ID), cell type, total biomass, colony area, Voronoi area, nearest neighbor, mean neighbor distance, etc.
    """
    if not hasattr(obj,'colonies'):
        obj.get_mothers()
    obj.count_colony_area(obj.Timesteps[-1])
    D_suc = obj.metadata['Diff_c']['suc']
    mu_ecw = obj.metadata['ecw']['GrowthRate']
    mu_cy = obj.metadata['cyano']['GrowthRate']
    df = obj.colonies.copy()
    # calculate voronoi area
    dfs = list()
    for type_ in df.type.unique():
        IDs = df[(df.Timestep ==0) & (df.type == type_)][['mother_cell','type']]
        points = df[(df.Timestep ==0) & (df.type == type_)][['x','y']].values
        vor = Voronoi(points)
        areas = [abs(np.sum( [0.5, -0.5] * vor.vertices[vor.regions[i]] * np.roll( np.roll(vor.vertices[vor.regions[i]], 1, axis=0), 1, axis=1) )) for i in range(len(vor.regions))]
        IDs.loc[:,'Voronoi Area'] = areas[1:]
        dfs.append(IDs)
    metrics = pd.concat(dfs)
    metrics.loc[:,'IPTG'] = obj.IPTG
    # total biomass per colony
    biomasses = df[df.Timestep==obj.Timesteps[-1]].groupby('mother_cell').sum().reset_index()[['mother_cell','biomass']]
    biomasses.columns=['mother_cell','total biomass']

    # Calculate nearest neighbors
    df3 = df[df.Timestep == 0]
    arr = df3[['x','y','z']].to_numpy()
    tree = KDTree(arr)
    d, i = tree.query(arr, k=2)
    arr1 = df3[df3.type==1][['x','y','z']].to_numpy()
    tree1 = KDTree(arr1)
    d1, i1 = tree1.query(df3[['x','y','z']].to_numpy(), k=2)
    arr2 = df3[df3.type==2][['x','y','z']].to_numpy()
    tree2 = KDTree(arr2)
    d2, i2 = tree2.query(df3[['x','y','z']].to_numpy(), k=2)
    n1 = list()
    n2 = list()
    for i in range(len(d1)):
        if d1[i,0]==0:
            n1.append(d1[i,1])
        else:
            n1.append(d1[i,0])
    for i in range(len(d2)):
        if d2[i,0]==0:
            n2.append(d2[i,1])
        else:
            n2.append(d2[i,0])
    df3.loc[:,'Nearest 1']=n1
    df3.loc[:,'Nearest 2']=n2
    df3.loc[:,'Nearest Neighbor'] = d[:,1]
    df3 = df3[['mother_cell','Nearest 1','Nearest 2','Nearest Neighbor']]
    df3.loc[:,'IC1'] = df3.loc[:,'Nearest 1'].mean()
    df3.loc[:,'IC2'] = df3.loc[:,'Nearest 2'].mean()
    df3.loc[:,'IC'] = df3.loc[:,'Nearest Neighbor'].mean()
    df3.loc[:,'Relative Neighbor Dist 1'] = df3.loc[:,'Nearest 1']/df3.loc[:,'IC1']
    df3.loc[:,'Relative Neighbor Dist 2'] = df3.loc[:,'Nearest 2']/df3.loc[:,'IC2']
    df3.loc[:,'Relative Neighbor Dist'] = df3.loc[:,'Nearest Neighbor']/df3.loc[:,'IC']
    df3.loc[:,'Z1'] = df3.loc[:,'Relative Neighbor Dist 1']/np.sqrt(D_suc/mu_cy)
    df3.loc[:,'Z2'] = df3.loc[:,'Relative Neighbor Dist 2']/np.sqrt(D_suc/mu_ecw)
    df3.loc[:,'Z1_2'] = df3.loc[:,'Relative Neighbor Dist 1']/np.sqrt(D_suc/mu_ecw)
    df3.loc[:,'Z2_1'] = df3.loc[:,'Relative Neighbor Dist 2']/np.sqrt(D_suc/mu_cy)
    df3.loc[:,'LogNearest 1'] = np.log(df3['Nearest 1'])
    df3.loc[:,'LogNearest 2'] = np.log(df3['Nearest 2'])
    df3.loc[:,'LogNearest'] = np.log(df3.loc[:,'Nearest Neighbor'])
    inv1 = list()
    for i in df[df.Timestep == 0][['x','y','z']].to_numpy():
        d1, i1 = tree1.query(i,k=2)
        if d1[0]==0:
            inv1.append(np.sum(1/d1[1]))
        else:
            inv1.append(np.sum(1/d1[0]))
    inv2 = list()
    inv3 = list()
    for i in df[df.Timestep == 0][['x','y','z']].to_numpy():
        d2, i2 = tree2.query(i,k=2)
        if d2[0]==0:
            inv2.append(np.sum(1/d2[1]))
        else:
            inv2.append(np.sum(1/d2[0]))
    # Calculate log inverse squared neighbor distance
    log_inv1 = list()
    for i in df[df.Timestep == 0][['x','y','z']].to_numpy():
        d1, i1 = tree1.query(i,k=2)
        if d1[0]==0:
            log_inv1.append(np.log(np.sum(1/(d1[1]**2))))
        else:
            log_inv1.append(np.log(np.sum(1/(d1[0]**2))))
    log_inv2 = list()
    for i in df[df.Timestep == 0][['x','y','z']].to_numpy():
        d2, i2 = tree2.query(i,k=2)
        if d2[0]==0:
            log_inv2.append(np.log(np.sum(1/(d2[1]**2))))
        else:
            log_inv2.append(np.log(np.sum(1/(d2[0]**2))))

    df3.loc[:,'Inv1']=inv1
    df3.loc[:,'Inv2']=inv2
    df3.loc[:,'Log Inv1']=log_inv1
    df3.loc[:,'Log Inv2']=log_inv2

    colony_area = df[df.Timestep==obj.Timesteps[-1]][['mother_cell','Colony Area']].drop_duplicates()

    metrics = metrics.merge(biomasses,on='mother_cell')
    metrics = metrics.merge(df3,on='mother_cell')
    metrics = metrics.merge(colony_area,on='mother_cell')
    metrics
    return metrics
