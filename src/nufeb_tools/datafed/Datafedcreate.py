import os
import getpass
import subprocess
from platform import platform
import sys
import argparse
import pickle
import json # For dealing with metadata
from datafed.CommandLib import API
import logging

from nufeb_tools import __version__

__author__ = "Jsakkos"
__copyright__ = "Jsakkos"
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
    parser.add_argument('--m', dest='metadata', action='store',
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

def main(args):
    """Warpper function

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--verbose", "42"]``).
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting crazy calculations...")
    filename = args.file
    file_title=args.name
    global_coll_id = args.id
    df_api = API()
    pkl_file = args.metadata

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
    _logger.info("Script ends here")


def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    """
    main(sys.argv[1:])


if __name__ == "__main__":

    run()

