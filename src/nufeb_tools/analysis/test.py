import os
import argparse
import sys
from matplotlib import pyplot as plt 
import h5py
import dask
import dask.array as da
import dask.dataframe as dd
from pathlib import Path
import pandas as pd
import seaborn as sns
p = Path('.')
folders = [x for x in p.iterdir() if x.is_dir()]


class NUFEB_data:
    """
    NUFEB simulation data class to collect results for analysis
    """
    def __init__(self):
        self.directory = r'D:\CADES Files\runs\Run_15_75_56_1'
        #self.id = args.id
        self.sucRatio = int(self.directory.split('_')[-2])
        if self.directory:
            self.get_local_data()
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
    def plot_overall_growth(self,ax=None, **kwargs):
        """
        This is a function to generate growth curve plots
        
        Args:
            ax (plt.ax)
            Axis to plot data on
            
            **kwargs
        """
        Biomass = self.biomass
        Biomass.index = Biomass.step/60/60*10 #convert timesteps (10s) to hours
        Biomass.index.name='Hours'
        Biomass.iloc[:,1:]=Biomass.iloc[:,1:]*1e18
        ax = ax or plt.gca()
        # Plot biomass over time
        ax.plot(Biomass.iloc[:,1],label='S. elongatus',color='#2ca25f')
        ax.plot(Biomass.iloc[:,2],label='E. coli', color ='#de2d26')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.set_axis_off()
        #axs[i].legend(frameon=False)
        ax.set_yscale('log')
        return ax
    def plot_average_nutrients(self,nutrient,ax=None,legend = None,**kwargs):
        sns.set_context('talk')
        sns.set_style('white')
        avgc = self.avg_con
        avgc.index = avgc.Time/60/60*10 #convert timesteps (10s) to hours
        avgc.index.name='Hours'
        avgc.drop('Time',inplace=True,axis=1)
        #avgc_f.iloc[:,1:]=Biomass.iloc[:,1:]*1e15
        SucroseMW = 342.3
        O2MW = 32
        CO2MW = 44.01
        avgc.O2 = avgc.O2/O2MW*1e3
        avgc.Sucrose = avgc.Sucrose/SucroseMW*1e3
        avgc['CO2'] = avgc['CO2']/CO2MW*1e3
        ax = ax or plt.gca()
        if nutrient == 'Sucrose':
            ax.plot(avgc.iloc[:,1],label='Sucrose',**kwargs)
        elif nutrient == 'CO2':
            ax.plot(avgc.iloc[:,2],label='CO2',**kwargs)
        elif nutrient == 'O2':
            ax.plot(avgc.iloc[:,0],label='O2',**kwargs)
        else:
            print('No nutrient specified')
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel('Concentration (mM)')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        if legend:
            ax.legend(frameon=False)
        return ax
    def radius_key(timestep):
        """
        generate the appropriate key for a radius at a given timestep
        Args:

        timestep (str)

        Example
        ts='1000'
        radius_key(ts)
        """
        return(f"radius{timestep}")
    


x = NUFEB_data()
f, axes = plt.subplots(ncols=3,nrows=2)
for ax in axes.ravel():
    x.plot_overall_growth(ax)
f, ax = plt.subplots()
sns.set_context('talk')
sns.set_style('white')
x.plot_average_nutrients('Sucrose',color='Green')

#%%
# 



print(numsteps)
hrs = [int(x)/360 for x in timepoints]
df = pd.DataFrame(columns=['id','type','time','biomass'])
# loop over all cell ids, c, from time = 0
for c in f['id']['0']:
    # loop over all timepoints t
    for t,h in zip(timepoints,hrs):
        # get list of ids from timepoint t
        ids = f['id'][t]
        arr = ids.__array__(ids.dtype)
        # make sure cell id matches id of interest
        i = np.where(arr == c)[0][0]
        # get radius
        radius = f[radius_key(t)][i]
        volume = 4/3*np.pi*radius**3 #volume in m^3
        # get celltype
        celltype=f['type'][t][i]
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

        #%% make growth curves
        rows=round(np.sqrt(len(df.id.unique())))
cols = int(np.ceil(len(df.id.unique())/rows))
fig ,axes = plt.subplots(nrows=rows,ncols=cols,sharex=True,figsize=(14,8))
axs = axes.ravel()
for i in df.id.unique():
    celltype=df[(df.id==i) & (df.time ==0)].type.values[0]
    if celltype==1:
        color = '#2ca25f'
    elif celltype ==2:
        color = '#de2d26'
    elif celltype ==0:
        print('Celltype is 0',i,c,celltype,t)
    axs[i-1].plot(df[(df.id==i)].biomass.values,c=color)
    #axs[i].set_title(f['type']['0'][c])
for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_axis_off()
    #ax.spines['bottom'].set_visible(False)
    
fig.tight_layout()
plt.show()
#%%
#plot cell size vs time
palette = sns.color_palette("mako_r", 6)
data = df[df.id==1].reset_index()
fig, ax = plt.subplots(figsize=(5,4))
ax.plot(data.time,data.biomass,color='#2ca25f')
#sns.lineplot(x='time',y='biomass',data=data,color = 'green',ax=ax)
for line in find_peaks(data.biomass):
    ax.vlines(data.time[line],data.biomass.min(),data.biomass.max()*1.1,color='#bdbdbd',ls=':')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Biomass (fg)')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
sns.despine()
fig.tight_layout()
fig.savefig('CyanogrowthHigh.png',dpi=600)
#plot cell size vs time
palette = sns.color_palette("mako_r", 6)
data = df[df.id==20].reset_index()
fig, ax = plt.subplots(figsize=(5,4))
ax.plot(data.time,data.biomass,color='#de2d26')
#sns.lineplot(x='time',y='biomass',data=data,color = 'green',ax=ax)
for line in find_peaks(data.biomass):
    ax.vlines(data.time[line],data.biomass.min(),data.biomass.max()*1.1,color='#bdbdbd',ls=':')
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Biomass (fg)')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
sns.despine()
fig.tight_layout()

#%% growth vs division
p0 = 0
for p1 in pks:
    #plt.plot(data.time[p0:p1],data.diameter[p0:p1])
    #plt.show()
    dt = data.time[p1]-data.time[p0]
    dy = data.biomass[p1]-data.biomass[p0]
    dydt = dy/dt
    #print(dydt)
    p0=p1+1

fig, axes = plt.subplots(ncols=2,figsize=(14,7))
celltypes = df.type.unique()
celltypes.sort()
for ct in celltypes:
    divs = pd.DataFrame(columns=['id','division','rate'])
    cells = df.id.unique()
    #cells.sort()
    for cell in cells:
        data = df[(df.id==cell) & (df.type==ct)].reset_index()
        pks,_ = find_peaks(data.biomass)
        p0 = 0
        for i,p1 in enumerate(pks):
            #plt.plot(data.time[p0:p1],data.diameter[p0:p1])
            #plt.show()
            dt = data.time[p1]-data.time[p0]
            dy = data.biomass[p1]-data.biomass[p0]
            dydt = dy/dt
            divs = divs.append(pd.DataFrame([[cell,i+1,dydt]],columns=['id','division','rate']),ignore_index=True)
            #print(dydt)
            p0=p1+1
    #plot cell id vs division rate over time
    piv = divs.pivot_table(index='id', columns='division', values='rate')
    g = sns.heatmap(piv, cmap='coolwarm',ax=axes[ct-1])
    cbar = g.collections[0].colorbar
    cbar.ax.set_ylabel(r'Growth rate ($\frac{fg}{hr}$)')
axes[0].set_title('S. elongatus')
axes[1].set_title('E. coli')
fig.tight_layout()

#%%
periods = 3
df['rate'] = df.diff(periods=periods).biomass/df.diff(periods=periods).time
df['rate'][df.rate < 0] = 0

sns.set_context('talk')
df.time = df.time.round(1)
fig, axes = plt.subplots(ncols=2,figsize=(14,7))
celltypes = df.type.unique()
celltypes.sort()
for ct in celltypes:
    rates = df[df.type==ct][['id','time','rate']]
    #plot cell id vs division rate over time
    piv = rates.pivot_table(index='id',columns='time', values='rate')
    g = sns.heatmap(piv, cmap='coolwarm',ax=axes[ct-1])
    cbar = g.collections[0].colorbar
    cbar.ax.set_ylabel(r'Growth rate ($\frac{fg}{hr}$)')
    #axes[ct-1].xaxis.set_major_locator(plt.MaxNLocator(5))
    axes[ct-1].set_xticklabels(axes[ct-1].get_xticklabels(), rotation=0) 
    axes[ct-1].locator_params(axis="x", nbins=6)
    axes[ct-1].set_xlabel('Time (hrs)')
    
axes[0].set_title('S. elongatus')
axes[1].set_title('E. coli')

fig.tight_layout()