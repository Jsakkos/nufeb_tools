from nufeb_tools import utils,plot
from scipy.integrate import odeint
import numpy as np
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import time
# import dask
# import dask.array as da
# import dask.dataframe as dd
from itertools import combinations
from scipy.spatial.distance import pdist,squareform

#%%
x = utils.get_data(test=True)
df = x.positions
#%%
timestep = 10000
combs = list(combinations(df[df.Timestep==timestep].ID,2))
df2 = df[df.Timestep == timestep].set_index(['ID'])
df2.sort_index(inplace=True)
# distances =pdist(df2[['x','y','z']])
pairwise = pd.DataFrame(
    squareform(pdist(df2[['x','y','z']])),
    columns = df2.index,
    index = df2.index
)
pairwise[pairwise ==0] = np.nan
vmin = pairwise.min().min()
sns.heatmap(pairwise,vmin=vmin)
# pdist(df[(df.Timestep==0) & (df.ID==1)][['x','y','z']].values,df[(df.Timestep==0) & (df.ID==2)][['x','y','z']].values)
#dist = np.sqrt(((df[(df.Timestep==0) & (df.ID==1)][['x','y','z']].squeeze() - df[(df.Timestep==0) & (df.ID==95)][['x','y','z']].squeeze())**2).sum())
# temp = (df[df.ID ==id][['x','y','z']].squeeze() - df[df.ID !=id][['x','y','z']])**2
# dist = pd.Series(np.sqrt(temp.x + temp.y + temp.z),name='Distance')
# pd.concat([df[df.ID !=id][['ID','type']],dist],axis=1).reset_index(drop=True)



# tic = time.perf_counter()

# dfs = list()
# for t in x.timepoints:
#     dfs.append(pd.concat([
#             pd.Series(np.ones(x.h5['x'][str(t)].len())*int(t),dtype=int,name='Timestep'),
#             pd.Series(x.h5['id'][str(t)],name='ID'),
#             pd.Series(x.h5['type'][str(t)],name='type'),
#             pd.Series(x.h5['x'][str(t)],name='x'),
#             pd.Series(x.h5['y'][str(t)],name='y'),
#             pd.Series(x.h5['z'][str(t)],name='z')],axis=1))
# df = pd.concat(dfs)
# toc = time.perf_counter()
# print(f"Ran in {toc - tic:0.4f} seconds")
#%%
# tic = time.perf_counter()
# dfs = list()
# for t in x.timepoints:
#     dfs.append(x.get_positions(t))
# df = pd.concat(dfs)
# toc = time.perf_counter()
# print(f"Ran in {toc - tic:0.4f} seconds")
#%% Dask version

# from dask import delayed
# tic = time.perf_counter()
# dfs = list()
# for t in x.timepoints:               # Use for loops to build up computation
#     dfs.append(delayed(x.get_positions)(t))     # Delay execution of function
# #dfs.compute()
# df = pd.concat(dfs)

# toc = time.perf_counter()