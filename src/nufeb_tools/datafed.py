import os
import subprocess
import sys
import argparse
import pickle
import json # For dealing with metadata
from datafed.CommandLib import API
import logging

from nufeb_tools import __version__

__author__ = "Jonathan Sakkos"
__copyright__ = "Jonathan Sakkos"
__license__ = "MIT"

_logger = logging.getLogger(__name__)
def parse_args(args):
    """Argument parser

    Args:
      args (List(str))
    """
    parser = argparse.ArgumentParser(description='Create datasets for DataFed')
    parser.add_argument('--id', dest='id', action='store',
                    help='Collection id to create data within')
    parser.add_argument('--f', dest='file', action='store',
                    help='Filename')
    parser.add_argument('--n', dest='name', action='store',
                    help='Name of dataset')
    parser.add_argument('--m', dest='metadata', action='store', default='',
                    help='Metadata file')
    args = parser.parse_args()
    parser.add_argument(
        "-v",
        "--verbose",
        dest="loglevel",
        help="set loglevel to INFO",
        action="store_const",
        const=logging.INFO,
    )
    parser.add_argument(
        "-vv",
        "--very-verbose",
        dest="loglevel",
        help="set loglevel to DEBUG",
        action="store_const",
        const=logging.DEBUG,
    )
    return parser.parse_args(args)


def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S"
    )

def upload(file, title, collection_id,metadata_file):
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
        _logger.info(pput_msg)
    else:
        _logger.debug('No metadata file found')
        sys.exit(1)

def create_collection(n_cyanos, n_ecw, SucPct,dims):
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
        df_api = API()
        df_api.setContext('p/eng107')
        collectionName = f'NUFEB_{n_cyanos}_{n_ecw}_{SucPct}_{dims[0]}_{dims[1]}_{dims[2]}'
        parent_collection = df_api.getAuthUser().split('/')[1]
        coll_msg = df_api.collectionCreate(collectionName,
                                        parent_id=parent_collection)
        global_coll_id = coll_msg[0].coll[0].id
        _logger.info(global_coll_id)
    except:
        global_coll_id = None
        _logger.debug('Unable to create collection')
    return global_coll_id

def verify_connection():
    """
    Verify Datafed installation
    """
    
    try:
        # This package is not part of anaconda and may need to be installed.
        from datafed.CommandLib import API
    except ImportError:
        _logger.info('datafed not found. Installing from pip.')
        subprocess.call([sys.executable, "-m", "pip", "install", 'datafed'])
        from datafed.CommandLib import API

    from datafed import version as df_ver

    if not df_ver.startswith('1.1.0'):
        _logger.info('Attempting to update DataFed.')
        subprocess.call([sys.executable, "-m", "pip", "install", '--upgrade', 'datafed'])
        _logger.debug('Please restart the python kernel or upgrade manually to V 1.1.0:1 if you are repeatedly seeing this message via'
            '\n\tpip install --upgrade datafed')
    else:
        df_api = API()
        #print('Success! You have DataFed: ' + df_ver)
    # Verify user authentication
    if not df_api.getAuthUser():
        _logger.info('You have not authenticated into DataFed Client')


    # Check default Globus endpoint
    if not df_api.endpointDefaultGet():
        endpoint = 'cades#CADES-OR'
        df_api.endpointDefaultSet(endpoint)



    #print('Your default Globus Endpoint in DataFed is:\n' + df_api.endpointDefaultGet())
    # Test the endpoint
    dget_resp = df_api.dataGet('d/35437908',
                            '/~/',
                            wait=True)
    _logger.debug(dget_resp)
    if dget_resp[0].task[0].status == 3:
        os.remove('35437908.md5sum')
        sys.exit(0)
    else:
        if dget_resp[0].task[0].msg == "globus connect offline":
            _logger.info('You need to activate your Globus Endpoint and/or ensure Globus Connect Personal is running.\n'
                'Please visit https://globus.org to activate your Endpoint')
            sys.exit(1)
        elif dget_resp[0].task[0].msg == "permission denied":
            _logger.info('Globus does not have write access to this directory. \n'
                'If you are using Globus Connect Personal, ensure that this notebook runs within'
                'one of the directories where Globus has write access. You may consider moving this'
                'notebook to a valid directory or add this directory to the Globus Connect Personal settings')
            sys.exit(1)
        else:
            NotImplementedError('Get in touch with us or consider looking online to find a solution to this problem:\n' + dget_resp[0].task[0].msg)
            sys.exit(1)
def main(args):
    """Wrapper function

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--verbose", "42"]``).
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Checking datafed connection")
    verify_connection()
 

def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    """
    main(sys.argv[1:])


if __name__ == "__main__":

    run()

