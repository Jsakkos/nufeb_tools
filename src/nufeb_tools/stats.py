
from scipy.spatial import KDTree
from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
import pandas as pd
from nufeb_tools import __version__

def fitness_metrics(obj):
    obj.count_colony_area(35000)
    df = obj.colonies.copy()
    # calculate voronoi area
    dfs = list()
    for type_ in df.type.unique():
        IDs = df[(df.Timestep ==0) & (df.type == type_)][['mother_cell','type']]
        points = df[(df.Timestep ==0) & (df.type == type_)][['x','y']].values
        vor = Voronoi(points)
        areas = [abs(np.sum( [0.5, -0.5] * vor.vertices[vor.regions[i]] * np.roll( np.roll(vor.vertices[vor.regions[i]], 1, axis=0), 1, axis=1) )) for i in range(len(vor.regions))]
        IDs['Voronoi Area'] = areas[1:]
        dfs.append(IDs)
    metrics = pd.concat(dfs)

    # total biomass per colony
    biomasses = df[df.Timestep==obj.Timesteps[-1]].groupby('mother_cell').sum().reset_index()[['mother_cell','biomass']]
    biomasses.columns=['mother_cell','total biomass']

    # Calculate nearest neighbors
    df3 = df[df.Timestep == 0]
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
    df3.loc[:,'Nearest1']=n1
    df3.loc[:,'Nearest2']=n2
    df3 = df3[['mother_cell','Nearest1','Nearest2']]
    # calculate sum of inverse neighbor distance
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

    #df[df.Timestep==0]
    metrics = metrics.merge(biomasses,on='mother_cell')#.groupby('mother_cell').max().reset_index()
    metrics = metrics.merge(df3,on='mother_cell')
    metrics = metrics.merge(colony_area,on='mother_cell')
    metrics
    return metrics
