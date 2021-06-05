import os
import argparse
import sys
from glob import glob
import h5py
import dask
import dask.array as da
import dask.dataframe as dd
from pathlib import Path

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
    def __init__(self, args):
        self.directory = args.directory
        self.id = args.id
        self.sucRatio = int(self.directory.split('_')[-2])
        if self.directory:
            self.get_local_data()
    def get_local_data(self):
        self.h5 = h5py.File(os.path.join(self.directory,'trajectory.h5'))
        self.biomass = dd.read_csv(os.path.join(self.directory,'Results','biomass.csv'))
        self.ntypes = dd.read_csv(os.path.join(self.directory,'Results','ntypes.csv'))
        self.avg_con = dd.read_csv(os.path.join(self.directory,'Results','avg_concentration.csv'))



def main(args):
    """Warpper function

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