import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def overall_growth(df,ax=None, **kwargs):
    """
    This is a function to generate growth curve plots
    
    Args:
        df (pd.DataFrame):
        Pandas Dataframe containing biomass data over time

        ax (plt.ax):
            Axis to plot data on
        
        **kwargs
    """
    Biomass = df
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
    
def average_nutrients(df,nutrient,ax=None,legend = None,**kwargs):
    """
    This is a function to plot the nutrient concentration over time
    
    Args:
        df (pd.DataFrame):
        Pandas Dataframe containing nutrient data

        nutrient (str):
            Name of the nutrient to plot, e.g., ``'Sucrose'``

        ax:
            Axis on which to make the plot

        legend (bool):
            Include legend in the plot

        **kwargs:
            Additional arguments to pass to plt.plot
    """
    sns.set_context('talk')
    sns.set_style('white')
    avgc = df
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