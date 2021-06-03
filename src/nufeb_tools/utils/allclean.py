"""
This is a script to remove old NUFEB simulation files and directories.
"""
import sys
import logging
import os
import shutil
import time
from glob import glob
import argparse
import logging
import sys
from pathlib import Path
from nufeb_tools import __version__

__author__ = "Jsakkos"
__copyright__ = "Jsakkos"
__license__ = "MIT"

_logger = logging.getLogger(__name__)

def parse_args(args):
    """Parse command line parameters

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--help"]``).

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(description="Delete old simulation files")
    parser.add_argument(
        "--version",
        action="version",
        version="nufeb_tools {ver}".format(ver=__version__),
    )
    parser.add_argument(
        '--dir', 
        dest="directory", 
        default=os.getcwd()
        help="Directory to look for simulation files", 
        type=str)
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


def clean_old_runs(args):
    """Delete files from old NUFEB simulations
    """
    
    files = [f for f_ in [Path(args.directory).rglob(e) for e in ('*.in', '*.h5','*.pkl')] for f in f_]
    for path in Path(args.directory).rglob('*.c'):
        print(path.name)
    #     rm runs/snapshot_*
    # rm runs/grid_*
    # rm runs/dump_*
    # rm runs/atom_*
    # rm runs/*.h5
    # rm runs/Inputscript_*
    # rm -rf runs/Results
    # rm -rf runs/Run*
    # rm runs/run*
    # rm runs/*.out
    # rm runs/output.lammps
    # rm runs/log.lammps
    # rm runs/*.pkl
    # rm runs/*.sh
    # rm runs/*.tar.gz
    # newpath = r'C:/Pytorch temp/' 
    # subdirs=['Labels','Data','Validation/Labels','Validation/Data','Weights']
    # if os.path.exists(newpath):
    #     shutil.rmtree(newpath)
    #     time.sleep(5)


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
    """Wrapper allowing :func:`fib` to be called with string arguments in a CLI fashion

    Instead of returning the value from :func:`fib`, it prints the result to the
    ``stdout`` in a nicely formatted message.

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``["--verbose", "42"]``).
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Starting crazy calculations...")
    print("The {}-th Fibonacci number is {}".format(args.n, fib(args.n)))
    _logger.info("Script ends here")


def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])


if __name__ == "__main__":
    # ^  This is a guard statement that will prevent the following code from
    #    being executed in the case someone imports this file instead of
    #    executing it as a script.
    #    https://docs.python.org/3/library/__main__.html

    # After installing your project with pip, users can also run your Python
    # modules as scripts via the ``-m`` flag, as defined in PEP 338::
    #
    #     python -m nufeb_tools.skeleton 42
    #
    run()