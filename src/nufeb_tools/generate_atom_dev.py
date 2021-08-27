"""
This is a script to seed NUFEB simulations
"""

import random
import argparse
import numpy as np
import pickle
from string import Template
import json # For dealing with metadata
import os # For file level operations
import time # For timing demonstrations
import datetime # To demonstrate conversion between date and time formats
from datafed.CommandLib import API
from itertools import combinations
import itertools
from collections import defaultdict
import sys
import logging
from importlib import resources
from pathlib import Path
from nufeb_tools import __version__

__author__ = "Jonathan Sakkos"
__copyright__ = "Jonathan Sakkos"
__license__ = "MIT"

_logger = logging.getLogger(__name__)

# TODO update lmp and atom templates
CellInfo = {'cyano': {'GrowthRate' : round(0.06/3600,7),
    'min_length' : 1e-6, 'max_length' : 5e-6, 'Diameter' : 1e-6, 'Density' : 370,
    'Inertia' : {'ixx' : 0, 'iyy' : 0, 'izz' : 9.2e-23, 'ixy' : 0, 'ixz' : 0, 'iyz' : 0},
        'K_s' : {'light' : 3.5e-4,'o2' : 2e-4, 'suc' : 1e-2,'co2' : 1.38e-4},
    'Yield' : 0.55,'Maintenance' : 0,'Decay' : 0},
        'ecw': {'GrowthRate' : 2.7e-04,
    'min_length' : 1.94e-6, 'max_length' : 2.72e-6, 'Diameter' : 0.73e-6,'Density' : 236,
    'Inertia' : {'ixx' : 0, 'iyy' : 0, 'izz' : 9.2e-23, 'ixy' : 0, 'ixz' : 0, 'iyz' : 0},
        'K_s' : {'light' : 0,'o2' : 1e-3, 'suc' : 3.6,'co2' : 5e-2},
    'Yield' : 0.43,'Maintenance' : 9.50e-7,'Decay' : 2e-5}
}   


# arguments to modify the conditions of the simulation seeding

def parse_args(args):
    """Parse command line parameters

    Args:
        args (List[str]): command line parameters as list of strings
            (for example  ``["--help"]``).

    Returns:
        :obj:`argparse.Namespace`: command line parameters namespace
    """
    parser = argparse.ArgumentParser(description='Create atom definition files')
    parser.add_argument('--n', dest='num', action='store',
                       default=1,
                       help='Create atom definition files for NUFEB with --n #files desired (default is 1)')
    parser.add_argument('--r', dest='reps', action='store',
                       default=1,
                       help='Number of replicates')
    parser.add_argument('--c',dest='culture_type',action='store',default='co',
                        help='Set culture conditions with --c (co-culture), --ax-c (cyano), --ax-e (e.coli)')
    parser.add_argument('--cells', dest='cells_init',action='store',default=50,
                        help='Number of total cells to initialize simulation with')
    parser.add_argument('--co2', dest='co2', action='store',
                       default=1e3,
                       help='Set initial CO2 concentration (mM)')
    parser.add_argument('--d', dest='dims', action='store', type=str,
                       default='1e-4,1e-4,1e-5',
                       help='Set simulation box dimensions (m)')
    parser.add_argument('--ts', dest='timestep', action='store',
                       default=600,
                       help='Timestep length')
    parser.add_argument('--t', dest='ntimesteps', action='store',
                       default=600,
                       help='Number of timesteps to run')
    parser.add_argument('--suc', dest='sucrose', action='store',
                       default=1e-19,
                       help='Set initial sucrose concentration (mM)')
    parser.add_argument('--grid', dest='grid', action='store',
                       default=2,
                       help='Diffusion grid density (um/grid)')
    parser.add_argument('--mono', dest='monolayer', action='store',
                      default=True,
                      help='Set seed generation to monolayer of cells')
    parser.add_argument('--u', dest='user', action='store',
                      help='CADES/CNMS user ID')
    parser.add_argument('--vtk',dest='vtk',action='store',default=True,
                        help='Enable VTK dump')
    parser.add_argument('--hdf5',dest='hdf',action='store',default=True,
                        help='Enable HDF5 dump')
    parser.add_argument('--img',dest='img',action='store',default=False,
                        help='Enable image dump')
    parser.add_argument('--mov',dest='movie',action='store',default=True,
                        help='Enable movie dump')
    parser.add_argument('--datafed', dest='datafed', action = 'store', default=False,
                        help='DataFed Upload')
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






class Nutrient:
    """
    Nutrient class to define the chemicals present in the simulation volume,
    their concentrations, and their properties
    """

    def __init__(self, c, d,mw,state,bc):
        self.concentration = c
        self.diffusion = d
        self.moleculuarWeight = mw
        self.state = state
        self.boundary = bc
        if self.moleculuarWeight is not None:
            self.concentrationNufeb = np.format_float_scientific(self.concentration*self.moleculuarWeight*1e-3,precision=1)
        else:
            self.concentrationNufeb = np.format_float_scientific(self.concentration,precision=1)

class Cell:
    """
    Bacteria object class
    """
    def __init__(self,Species,Group,idx, args,Info = CellInfo):
        self.group = Group
        self.species = Species
        self.growth = Info[self.species]['GrowthRate']
        self.min_length = Info[self.species]['min_length']
        self.max_length = Info[self.species]['max_length']
        self.diameter = Info[self.species]['Diameter']
        self.density = Info[self.species]['Density']
        self.Ks = Info[self.species]['K_s']
        self.yld = Info[self.species]['Yield']
        self.maintenance = Info[self.species]['Maintenance']
        self.decay = Info[self.species]['Decay']
        self.inertia = Info[self.species]['Inertia']
        self.boundaries = [float(x) for x in args.dims.split(',')]
        self.length = random.uniform(self.min_length,self.max_length)
        self.monolayer = args.monolayer
        self.x = random.uniform(0,self.boundaries[0])
        self.y = random.uniform(0,self.boundaries[1])
        if self.monolayer:
            self.z = 1e-6
        else:
            self.z = random.uniform(0,self.boundaries[2])
        self.angle = random.uniform(1,180)
        self.x_displacement = np.format_float_scientific(self.length/2*np.cos(self.angle*np.pi/360),precision=1)
        self.y_displacement = np.format_float_scientific(self.length/2*np.sin(self.angle*np.pi/360),precision=1)
        self.z_displacement = 0
        self.z_angle = 0
        self.index = idx

    def Atom(self):
        """
        Function to return atom (cell) positions to render atom.in file
        """
        return ' '.join(map(str, [self.index+1,self.group, 1, self.density,
                np.format_float_scientific(self.x,precision=2),
                np.format_float_scientific(self.y,precision=2),
                np.format_float_scientific(self.z,precision=2),1, ' \n']))


    def rotate(self,z_dim = False):
        """
        Randomly generate cell orientation displacements based on input angle
        """
        self.x_displacement = np.format_float_scientific(self.length/2*np.cos(self.angle*np.pi/360),precision=1)
        self.y_displacement = np.format_float_scientific(self.length/2*np.sin(self.angle*np.pi/360),precision=1)
        zd = z_dim
        if zd == True:
            self.z_angle = random.uniform(1,180)
            self.z_displacement = np.format_float_scientific(self.length/2*np.sin(self.z_angle*np.pi/360),precision=1)
        else:
            self.z_displacement = 0
            self.z_angle = 0
        return [self.x_displacement, self.y_displacement, self.z_displacement]

    def Bacillus(self):
        """
        Function to return rod shape (bacillus) parameters to render atom.in file
        """
        return ' '.join(map(str, [self.index+1] + list(self.inertia.values()) + self.rotate() + [self.diameter] + [' \n']))

    def Report(self):
        """
        Return cell position, orientation, and size
        """
        return [self.x,self.y,self.z,self.angle,self.length,self.diameter]
    def Check(self):
        """
        Return cell position, orientation, and size
        """
        return self.x,self.y,self.length, self.index


class Culture:
    """
    Create a collection of cells with defined positions, lengths, and orientations
    """
    def __init__(self, args):
        self.cultureType = args.culture_type
        self.n_cells = int(args.cells_init)
        self.SucRatio = round(random.random(),3)
        self.SucPct = int(self.SucRatio*100)
        self.boundaries = [float(x) for x in args.dims.split(',')]
        if self.cultureType == 'co':
            self.cell_types = ['cyano','ecw']
            self.n_cyanos = int(random.uniform(1,self.n_cells-1))
            self.n_ecw = int(self.n_cells-self.n_cyanos)
            # self.n_cells = n_cyanos + n_ecw
            self.cellCount = {'cyano' : self.n_cyanos,'ecw' : self.n_ecw}
            self.cyGroup = 'group CYANO type 1'
            self.ecwGroup = 'group ECW type 2'
            self.cyDiv = f'fix div_cya CYANO nufeb/divide/bacillus {CellInfo["cyano"]["max_length"]} {random.randint(1,1e6)}'
            self.ecwDiv = f'fix div_ecw ECW nufeb/divide/bacillus {CellInfo["ecw"]["max_length"]} {random.randint(1,1e6)}'
            self.cyMonod = f'fix monod_cyano CYANO nufeb/monod/cyano light {CellInfo["cyano"]["K_s"]["light"] : .2e} o2 co2 {CellInfo["cyano"]["K_s"]["co2"] : .2e} suc gco2 growth {CellInfo["cyano"]["GrowthRate"] : .2e} yield {CellInfo["cyano"]["Yield"] : .2e} suc_exp {self.SucRatio}'
            self.ecwMonod = f'fix monod_ecw ECW nufeb/monod/ecoli/wild suc {CellInfo["ecw"]["K_s"]["suc"] : .2e} o2 {CellInfo["ecw"]["K_s"]["o2"] : .2e} co2 growth {CellInfo["ecw"]["GrowthRate"] : .2e} yield {CellInfo["ecw"]["Yield"] : .2e} maintain {CellInfo["ecw"]["Maintenance"] : .2e} decay {CellInfo["ecw"]["Decay"] : .2e}'
            self.cyanoCount = 'variable ncyano equal "count(CYANO)"'
            self.ecwCount = 'variable necw equal "count(ECW)"'
            self.vcyano = 'v_ncyano'
            self.vecw = 'v_necw'
        elif self.cultureType == 'ax-c':
            self.cell_types = ['cyano']
            self.n_cyanos = int(random.uniform(1,self.n_cells))
            self.n_ecw = 0
            # self.n_cells = n_cyanos
            self.cellCount = {'cyano' : self.n_cyanos}
            self.cyGroup = 'group CYANO type 1'
            self.ecwGroup = ''
            self.cyDiv = f'fix div_cya CYANO nufeb/divide/bacillus {CellInfo["cyano"]["max_length"]} {random.randint(1,1e6)}'
            self.ecwDiv = ''
            self.cyMonod = f'fix monod_cyano CYANO nufeb/monod/cyano light {CellInfo["cyano"]["K_s"]["light"] : .2e} o2 co2 {CellInfo["cyano"]["K_s"]["co2"] : .2e} suc gco2 growth {CellInfo["cyano"]["GrowthRate"] : .2e} yield {CellInfo["cyano"]["Yield"] : .2e} suc_exp {self.SucRatio}'
            self.ecwMonod = ''
            self.cyanoCount = 'variable ncyano equal "count(CYANO)"'
            self.ecwCount = ''
            self.vcyano = 'v_ncyano'
            self.vecw = ''
        elif self.cultureType == 'ax-e':
            self.cell_types = ['ecw']
            self.n_ecw = int(random.uniform(1,self.n_cells))
            self.n_cyanos=0
            # self.n_cells = n_ecw
            self.cellCount = {'ecw' : self.n_ecw}
            self.cyGroup = ''
            self.ecwGroup = 'group ECW type 1'
            self.cyDiv = ''
            self.ecwDiv = f'fix div_ecw ECW nufeb/divide/bacillus {CellInfo["ecw"]["max_length"]} {random.randint(1,1e6)}'
            self.cyMonod = ''
            self.ecwMonod = f'fix monod_ecw ECW nufeb/monod/ecoli/wild suc {CellInfo["ecw"]["K_s"]["suc"] : .2e} o2 {CellInfo["ecw"]["K_s"]["o2"] : .2e} co2 growth {CellInfo["ecw"]["GrowthRate"] : .2e} yield {CellInfo["ecw"]["Yield"] : .2e} maintain {CellInfo["ecw"]["Maintenance"] : .2e} decay {CellInfo["ecw"]["decay"] : .2e}'
            self.cyanoCount = ''
            self.ecwCount = 'variable necw equal "count(ECW)"'
            self.vcyano = ''
            self.vecw = 'v_necw'
        counter = itertools.count(0)
        self.cells = [Cell(CellType,c,next(counter),args) for c, CellType in enumerate(self.cell_types,start=1) for i in range(self.cellCount[CellType])]
        self.Check()
    def __iter__(self):
        return iter(self.cells)

    def Check(self):
        """
        Check initial positions of each cell and move one of them if there is a collision
        """
        cleared = False
        while not cleared:
            for i in list(combinations([cell.Check() for cell in self.cells],2)):
            # for i in list(combinations(zip(self.locations.x,self.locations.y,self.locations.length,self.locations.index),2)):
                x1 = i[0][0]
                y1 = i[0][1]
                r1 = i[0][2]/2
                idx1 = i[0][3]
                x2 = i[1][0]
                y2 = i[1][1]
                r2 = i[1][2]/2
                idx1 = i[0][3]
                idx2 = i[1][3]
                distance = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2);
                radii = (r1 + r2) * (r1 + r2);
                if (distance == radii):
                    cleared = True
                elif (distance > radii):
                    cleared = True
                else:
                    if x1 > x2 and y1 > y2:
                        if x1 + r1 > 0 and x1 + r1 < self.boundaries[0] and y1 + r1 > 0 and y1 + r1 < self.boundaries[1]:
                            self.cells[idx1].x = x1 + r1/2
                            self.cells[idx1].y = y1 + r1/2
                    elif x1 > x2 and y1 < y2:
                        if x1 + r1 > 0 and x1 + r1 < self.boundaries[0] and y1 - r1 > 0 and y1 - r1 < self.boundaries[1]:
                            self.cells[idx1].x = x1 + r1/2
                            self.cells[idx1].y = y1 - r1/2
                    elif x1 < x2 and y1 > y2:
                        if x1 - r1 > 0 and x1 - r1 < self.boundaries[0] and y1 + r1 > 0 and y1 + r1 < self.boundaries[1]:
                            self.cells[idx1].x = x1 - r1/2
                            self.cells[idx1].y = y1 + r1/2
                    else:
                        if x1 - r1 > 0 and x1 - r1 < self.boundaries[0] and y1 - r1 > 0 and y1 - r1 < self.boundaries[1]:
                            self.cells[idx1].x = x1 - r1/2
                            self.cells[idx1].y = y1 - r1/2
                        _logger.debug(f'Bumped from {x1 :.2e}, {y1 :.2e} to {self.cells[idx1].x :.2e}, {self.cells[idx1].y :.2e}')
                    cleared = False
        return

def clean():
    if os.path.isdir('runs'):
        import shutil
        try:
            shutil.rmtree('runs')
        except OSError as e:
            print("Error: %s : %s" % ('runs', e.strerror))


def main(args):
    """Wrapper function to generate new NUFEB simulation conditions

    Args:
      args (List[str]): command line parameters as list of strings
          
    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.debug("Generating NUFEB simulation files")
    
    # create nutrients
    light = Nutrient(1e-1,None,None,'g','nn')
    co2 = Nutrient(float(args.co2),1.9e-09,44.01,'l','nn')
    o2 = Nutrient(0.28125,2.30e-9,32,'l','nn')
    sucrose = Nutrient(float(args.sucrose),5.2e-10,342.3,'l','nn')
    gco2 = Nutrient(0,None,44.01,'g','nn')
    TEMPLATES_DIR = (Path(__file__).parent) / 'templates'

    captureRate = round(1000/args.timestep)
    # define dump parameters
    dump_list = {'vtk_dump' : f'dump atom_vtk all vtk {captureRate} dump*.vtu id type diameter vx vy vz fx fy fz  \n dump grid_vtk all grid/vtk {captureRate} dump_%_*.vti con',
            'image_dump' : f'dump du_image all image {captureRate} image.*.jpg type diameter zoom 2 bacillus type size 1280 720 view 45 60 \n dump_modify du_image acolor 1 green acolor 2 red',
            'movie_dump' : f'dump du_mov all movie {captureRate} movie.avi type diameter zoom 1.5 bacillus type size 1280 720 view 0 0 \n dump_modify du_mov acolor 1 green acolor 2 red',
            'hdf_dump' : f'dump du_h5 all nufeb/hdf5 {captureRate} dump.h5 id type x y z vx vy vz fx fy fz radius conc reac'
            }


    dumps = defaultdict(list)
    for i in range(4):
        tmp = ['vtk_dump','image_dump','movie_dump','hdf_dump']
        dumps[tmp[i]]

    for dump,dump_var in zip([args.vtk,args.img,args.movie,args.hdf],['vtk_dump','image_dump','movie_dump','hdf_dump']):
        if dump is True or dump == 'True':
            dumps[dump_var] = dump_list[dump_var]
        else:
            dumps[dump_var] = ''

    ## Species-specific parameters



    # check for runs folder
    if not os.path.isdir('runs'):
        os.mkdir('runs')
    x = float(args.dims.split(',')[0])
    y = float(args.dims.split(',')[1])
    z = float(args.dims.split(',')[2])
    for n in range(1,int(args.num)+1):
        culture = Culture(args)
        atoms_list = []
        bacilli_list = []
        # Create list of atoms and bacilli for atom definition file
        for cell in culture.cells:
            atoms_list.append(cell.Atom())
            bacilli_list.append(cell.Bacillus())
        # make atom definition file
        for r in range(1,int(args.reps)+1):
            L = [' NUFEB Simulation\r\n\n',f'     {args.cells_init} atoms \n',
                f'     {len(culture.cell_types)} atom types \n',f'     {args.cells_init} bacilli \n\n',
                f'  0.0e-4   {x :.2e}  xlo xhi \n',f'  0.0e-4   {y :.2e}  ylo yhi \n',
                f'  0.0e-4   {z :.2e}  zlo zhi \n\n', ' Atoms \n\n'
                ]
            atoms = L+atoms_list
            atoms.append('\n')
            atoms.append(' Bacilli \n\n')
            atoms = atoms + bacilli_list
            #write atom definition file
            f= open(f"runs/atom_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}_{r}.in","w+")
            f.writelines(atoms)
        RUN_DIR = Path('runs') / f'Run_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}_{args.reps}'
        if not os.path.isdir(RUN_DIR):
            os.mkdir(RUN_DIR)
            #os.mkdir(f'runs/Run_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}_{args.reps}')
        #write initial conditions json file
        dumpfile = open(RUN_DIR / 'metadata.json','w')
        #dumpfile = open(f"/runs/Run_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}_{args.reps}/metadata.json",'w')
        json.dump(CellInfo, dumpfile, indent = 6)
        dumpfile.close()
        ###

        #write Inputscript
        #open the file
        filein = open( TEMPLATES_DIR / 'bacillus.txt' )
        #filein = resources.read_text("nufeb_tools.templates", "Bacillus.txt")
        #read it
        src = Template( filein.read() )
        #do the substitution
        result = src.safe_substitute({'n' : args.cells_init,
                                    'SucRatio' : culture.SucRatio,
                                    'SucPct' : culture.SucPct,
                                    'n_cyanos' : culture.n_cyanos,
                                    'n_ecw' : culture.n_ecw,
                                    'Replicates' : args.reps,
                                    'Timesteps' : args.ntimesteps,
                                    'ts' : args.timestep,
                                    'CYANOGroup' : culture.cyGroup,
                                    'ECWGroup' : culture.ecwGroup,
                                    'Zheight' : float(args.dims.split(',')[2]),
                                    'CYANODiv'  : culture.cyDiv,
                                    'ECWDiv' : culture.ecwDiv,
                                    'light' : light.concentration,
                                    'co2' : co2.concentration,
                                    'o2' : o2.concentration,
                                    'sucrose' : sucrose.concentration,
                                    'gco2' : gco2.concentration,
                                    'CYANOMonod' : culture.cyMonod,
                                    'ECWMonod' : culture.ecwMonod,
                                    'CYANOcount' : culture.cyanoCount,
                                    'ECWcount' : culture.ecwCount,
                                    'v_ncyano' : culture.vcyano,
                                    'v_necw' : culture.vecw,
                                    'vtk_dump': dumps['vtk_dump'],
                                    'image_dump' : dumps['image_dump'],
                                    'movie_dump' : dumps['movie_dump'],
                                    'hdf_dump' : dumps['hdf_dump']
                                    })
        f= open(f"./runs/Inputscript_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}.lmp","w+")
        f.writelines(result)



        if args.datafed is True or args.datafed == 'True':
        #create DataFed collection to hold the results
            df_api = API()
            df_api.setContext('p/eng107')
            collectionName = f'NUFEB_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}_{x}_{y}_{z}'
            parent_collection = df_api.getAuthUser().split('/')[1]
            coll_msg = df_api.collectionCreate(collectionName,
                                            parent_id=parent_collection)
            global_coll_id = coll_msg[0].coll[0].id
        else:
            global_coll_id = None


        #write local run script
        #open the file
        filein = open( TEMPLATES_DIR / 'local.txt' )
        #filein = resources.read_text("nufeb_tools.templates", "local.txt")
        #read it
        src = Template( filein.read() )
        #do the substitution
        result = src.safe_substitute({'n' : n, 'SucRatio' : culture.SucRatio, 'SucPct' : culture.SucPct,
                                    'n_cyanos' : culture.n_cyanos, 'n_ecw' : culture.n_ecw,
                                    'Reps' : args.reps,'id': global_coll_id})
        f= open(f"./runs/local_{culture.n_cyanos}_{culture.n_ecw}_{culture.SucPct}.sh","w+")
        f.writelines(result)
    #write slurm script
    #open the file
    filein = open( TEMPLATES_DIR / 'slurm_dev.txt' )
    #filein = resources.read_text("nufeb_tools.templates", "Slurm.txt")
    #read it
    src = Template( filein.read() )
    #do the substitution
    result = src.safe_substitute({'n' : args.cells_init, 'job' : f"NUFEB_cyano{n}",
                                    'USER' : args.user,'Replicates'  : args.reps,
                                    'SucPct' : culture.SucPct,'n_cyanos' : culture.n_cyanos,
                                    'n_ecw' : culture.n_ecw,'id': global_coll_id})
    # TODO update slurm script
    #f= open(f"./RunBatch.slurm","w+")
    #f.writelines(result)
    _logger.info("Script ends here")

def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])


if __name__ == "__main__":

    run()
