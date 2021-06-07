import os
import argparse
import sys
from glob import glob
import h5py
#import dask
#import dask.array as da
#import dask.dataframe as dd
from pathlib import Path
import seaborn as sns
import pandas as pd
simulation_list = next(os.walk('D:\\CADES Files\\runs\\'))[1]
def parse_args(args):
    """Argument parser

    Args:
      args (List(str))
    """
    parser = argparse.ArgumentParser(description='Get datasets')
    parser.add_argument('--id', dest='id', action='store',
                    help='Collection id to get data from')
    parser.add_argument('--dir', dest='directory', action='store',
                    help='Local directory to look for simulation data')

    return parser.parse_args(args)

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
            ax:
                Axis to plot data on

            **kwargs:
                Additional arguments to pass to plt.plot
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
        """
        Function to plot the average nutrient concentration in the simulation volume over time

        Args:
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
        Generate the appropriate key for a radius at a given timestep

        Args:
            timestep (str):
                The simulation timestep to get keys from

        Example:

        >>> ts='1000'
        >>> radius_key(ts)
        """
        return(f"radius{timestep}")



def main(args):
    """Wrapper function

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--verbose", "42"]``).
    """
    args = parse_args(args)

    return NUFEB_data(args)
    
def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    """
    main(sys.argv[1:])


if __name__ == "__main__":

    run()