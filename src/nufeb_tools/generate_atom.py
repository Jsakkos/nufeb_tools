import random
from scipy.stats import loguniform
import argparse
import numpy as np
from string import Template
import json
import os
import sys
from datetime import date
import logging
from nufeb_tools import __version__
from pathlib import Path
from glob import glob


__author__ = "Jonathan Sakkos"
__copyright__ = "Jonathan Sakkos"
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
    # arguments to modify the conditions of the simulation seeding
    parser = argparse.ArgumentParser(description="Create atom definition files")
    parser.add_argument(
        "--n",
        dest="num",
        action="store",
        default=1,
        help="Create atom definition files for NUFEB with --n #files desired (default is 1)",
    )
    parser.add_argument(
        "--r", dest="reps", action="store", default=1, help="Number of replicates"
    )
    parser.add_argument(
        "--c",
        dest="culture_type",
        action="store",
        default="co",
        help="Set culture conditions with --c (co-culture), --ax-c (cyano), --ax-e (e.coli)",
    )
    parser.add_argument(
        "--co2",
        dest="co2",
        action="store",
        default=6.8e-1,
        help="Set initial CO2 concentration (mM)",
    )
    parser.add_argument(
        "--d",
        dest="dims",
        action="store",
        type=str,
        default="1e-4,1e-4,1e-5",
        help="Set simulation box dimensions (m)",
    )
    parser.add_argument(
        "--t",
        dest="timesteps",
        action="store",
        default=35000,
        help="Number of timesteps to run",
    )
    parser.add_argument(
        "--suc",
        dest="sucrose",
        action="store",
        default=1e-19,
        help="Set initial sucrose concentration (mM)",
    )
    parser.add_argument(
        "--grid",
        dest="grid",
        action="store",
        default=2,
        help="Diffusion grid density (um/grid)",
    )
    parser.add_argument(
        "--mono",
        dest="monolayer",
        action="store",
        default=True,
        help="Set seed generation to monolayer of cells",
    )
    parser.add_argument("--u", dest="user", action="store", help="CADES/CNMS user ID")
    parser.add_argument(
        "--cells",
        dest="cells",
        action="store",
        default=None,
        help="Number of cyanobacteria and e.coli to initialize simulation with, `e.g., 100,100. ` Default is random number between 1 and 100.",
    )
    parser.add_argument(
        "--od",
        dest="od",
        action="store",
        default=None,
        help="Optical density of cyanobacteria and e.coli to initialize simulation with, `e.g., 0.3,1`",
    )
    parser.add_argument(
        "--sucR",
        dest="SucRatio",
        action="store",
        default=None,
        help="Set sucrose secretion ratio (0 to 1). Default is random.",
    )
    parser.add_argument(
        "--iptg",
        dest="iptg",
        action="store",
        default=None,
        type=float,
        help="Set IPTG induction for sucrose secretion (0 to 1). Default is random.",
    )
    parser.add_argument(
        "--muecw",
        dest="mu_ecw",
        action="store",
        default=3.55e-4,
        type=float,
        help="E. coli W maximum growth rate",
    )
    parser.add_argument(
        "--mucya",
        dest="mu_cya",
        action="store",
        default=1.802e-5,
        type=float,
        help="S. elongatus maximum growth rate",
    )
    parser.add_argument(
        "--rhoecw",
        dest="rho_ecw",
        action="store",
        default=230,
        type=float,
        help="E. coli W cell density",
    )
    parser.add_argument(
        "--rhocya",
        dest="rho_cya",
        action="store",
        default=370,
        type=float,
        help="S. elongatus cell density",
    )
    parser.add_argument(
        "--kco2",
        dest="kco2",
        action="store",
        default=0.00812672,
        type=float,
        help="S. elongatus Kco2",
    )
    parser.add_argument(
        "--ksuc",
        dest="ksuc",
        action="store",
        default=3.6,
        type=float,
        help="E. coli W Ksuc",
    )
    parser.add_argument(
        "--maintecw",
        dest="ecw_maint",
        action="store",
        default=0,
        type=float,
        help="E. coli W maintenance cost",
    )
    parser.add_argument(
        "--yieldecw",
        dest="ecw_yield",
        action="store",
        default=0.49,
        type=float,
        help="E. coli W biomass yield (g-dw/g-sucrose)",
    )
    parser.add_argument(
        "--decayecw",
        dest="ecw_decay",
        action="store",
        default=0,
        type=float,
        help="E. coli W decay rate",
    )
    parser.add_argument(
        "--biodt",
        dest="biodt",
        action="store",
        default=100,
        type=int,
        help="Biological timesteps.",
    )
    parser.add_argument(
        "--mass",
        dest="mass_max",
        action="store",
        default=None,
        type=float,
        help="Maximum biomass",
    )

    parser.add_argument(
        "--vtk",
        dest="vtk",
        action="store",
        default=False,
        type=bool,
        help="Output VTK files",
    )
    parser.add_argument(
        "--h5", dest="hdf5", action="store", default=True, help="Output HDF5 files"
    )
    parser.add_argument(
        "--lammps",
        dest="lammps",
        action="store",
        default=False,
        help="Output lammps files",
    )
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
    parser.add_argument(
        "--niter",
        dest="niter",
        action="store",
        help="Number of iterations for diffusion calculation",
        type=int,
        default=1000000,
    )
    parser.add_argument(
        "--suc_halt",
        dest="sucrose_halt",
        default=None,
        type=int,
        help="Add halt condition for sucrose levels.",
    )
    parser.add_argument(
        "--division", dest="division", default="on", type=str, help="Cellular division"
    )
    return parser.parse_args(args)


# TODO Change sucRatio to IPTG
def setup_logging(loglevel):
    """Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    """
    logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s"
    logging.basicConfig(
        level=loglevel, stream=sys.stdout, format=logformat, datefmt="%Y-%m-%d %H:%M:%S"
    )


def clean():
    """Remove old NUFEB runs"""
    if os.path.isdir("runs"):
        import shutil

        # TODO add folder removal to logger
        try:
            shutil.rmtree("runs")
        except OSError as e:
            print("Error: %s : %s" % ("runs", e.strerror))
    slurm_path = glob("*.slurm")
    if slurm_path:
        for file in slurm_path:
            os.remove(file)


def main(args):
    """Wrapper function to generate new NUFEB simulation conditions

    Args:
      args (List[str]): command line parameters as list of strings

    """
    args = parse_args(args)
    setup_logging(args.loglevel)
    _logger.info("Generating NUFEB simulation files")
    # maximum growth rates, mu
    # $mu_cyanos = round(0.06/3600,7)
    # mu_ecw = 2.7e-04 for 37C only
    # mu_ecw = 6.71e-5
    # molecular weights of co2 and sucrose for unit conversions
    CO2MW = 44.01
    SucMW = 342.3
    Biomass2OD = np.array([0.28, 0.44]) * np.prod(
        [float(x) for x in args.dims.split(",")]
    )  # kg-cells/m^3/OD
    TEMPLATES_DIR = (Path(__file__).parent) / "templates"
    InitialConditions = {
        "cyano": {
            "GrowthRate": args.mu_cya,
            "min_size": 1.37e-6,
            "max_size": 1.94e-6,
            "Density": args.rho_cya,
            "K_s": {"sub": 3.5e-4, "o2": 2e-4, "suc": 1e-2, "co2": args.kco2},
            "GrowthParams": {"Yield": 0.55, "Maintenance": 0, "Decay": 0},
        },
        "ecw": {
            "GrowthRate": args.mu_ecw,
            "min_size": 8.8e-7,
            "max_size": 1.39e-6,
            "Density": args.rho_ecw,
            "K_s": {"sub": 0, "o2": 1e-3, "suc": args.ksuc, "co2": 5e-2},
            "GrowthParams": {
                "Yield": args.ecw_yield,
                "Maintenance": args.ecw_maint,
                "Decay": args.ecw_decay,
            },
        },
        "Nutrients": {
            "Concentration": {
                "sub": 1e-1,
                "o2": 9e-3,
                "suc": float(args.sucrose) * SucMW * 1e-3,
                "co2": float(args.co2) * CO2MW * 1e-3,
            },
            "State": {"sub": "g", "o2": "l", "suc": "l", "co2": "l"},
            "xbc": {"sub": "nn", "o2": "nn", "suc": "nn", "co2": "nn"},
            "ybc": {"sub": "nn", "o2": "nn", "suc": "nn", "co2": "nn"},
            "zbc": {"sub": "nn", "o2": "nd", "suc": "nn", "co2": "nd"},
        },
        "Diff_c": {"sub": 0, "o2": 2.30e-9, "suc": 5.2e-10, "co2": 1.9e-09},
        "Dimensions": [float(x) for x in args.dims.split(",")],
        "Replicates": int(args.reps),
    }
    # check for runs folder
    if not os.path.isdir("runs"):
        os.mkdir("runs")
    today = str(date.today())
    if args.mass_max is not None:
        max_mass = float(args.mass_max)
    else:
        max_mass = (
            1.5e-11
            * np.prod([float(x) for x in args.dims.split(",")])
            / (1e-4 * 1e-4 * 1e-5)
        )
    # TODO add file generation details to logger
    for n in range(1, int(args.num) + 1):
        if args.iptg is not None:
            IPTG = float(args.iptg)
            SucRatio = IPTG

        else:
            IPTG = np.round(loguniform.rvs(1e-3, 1e0, size=1)[0], 5)
            SucRatio = IPTG
        InitialConditions["IPTG"] = IPTG
        InitialConditions["SucRatio"] = SucRatio
        _logger.info(f"IPTG = {IPTG}")
        SucPct = int(SucRatio * 100)
        mean_cyano_mass = (
            np.mean(
                [
                    (InitialConditions["cyano"]["min_size"]) ** 3,
                    (InitialConditions["cyano"]["max_size"]) ** 3,
                ]
            )
            / 6
            * np.pi
            * 370
        )
        mean_ecw_mass = (
            np.mean(
                [
                    (InitialConditions["ecw"]["min_size"]) ** 3,
                    (InitialConditions["ecw"]["max_size"]) ** 3,
                ]
            )
            / 6
            * np.pi
            * 370
        )
        if args.culture_type == "co":
            cell_types = ["cyano", "ecw"]
            if args.cells is not None:
                n_cyanos = int(args.cells.split(",")[0])
                n_ecw = int(args.cells.split(",")[1])
                total_cyano_biomass = n_cyanos * mean_cyano_mass
                total_ecw_biomass = n_ecw * mean_ecw_mass
                InitialConditions["cyano"]["OD"] = total_cyano_biomass / Biomass2OD[0]
                InitialConditions["ecw"]["OD"] = total_ecw_biomass / Biomass2OD[1]
            elif args.od is not None:
                total_cyano_biomass = float(args.od.split(",")[0]) * Biomass2OD[0]
                n_cyanos = round(total_cyano_biomass / mean_cyano_mass)
                total_ecw_biomass = float(args.od.split(",")[1]) * Biomass2OD[1]
                n_ecw = round(total_ecw_biomass / mean_ecw_mass)
                InitialConditions["cyano"]["OD"] = float(args.od.split(",")[0])
                InitialConditions["ecw"]["OD"] = float(args.od.split(",")[1])

            else:
                n_cyanos = int(random.uniform(1, 100))
                n_ecw = int(random.uniform(1, 100))
                total_cyano_biomass = n_cyanos * mean_cyano_mass
                total_ecw_biomass = n_ecw * mean_ecw_mass
                InitialConditions["cyano"]["OD"] = total_cyano_biomass / Biomass2OD[0]
                InitialConditions["ecw"]["OD"] = total_ecw_biomass / Biomass2OD[1]
            n_cells = n_cyanos + n_ecw
            cyGroup = "group CYANO type 1"
            ecwGroup = "group ECW type 2"
            if args.division.lower() == "on":
                cyDiv = f"fix d1 CYANO divide {args.biodt} v_EPSdens v_divDia1 {random.randint(1,1e6)}"
                ecwDiv = f"fix d2 ECW divide {args.biodt} v_EPSdens v_divDia2 {random.randint(1,1e6)}"
            else:
                cyDiv = ""
                ecwDiv = ""
            masses = "c_myMass[1]+c_myMass[2]"

        elif args.culture_type == "ax-c":
            cell_types = ["cyano"]
            if args.cells is not None:
                n_cyanos = int(args.cells.split(",")[0])
            elif args.od is not None:
                total_cyano_biomass = float(args.od.split(",")[0]) * Biomass2OD[0]
                n_cyanos = round(total_cyano_biomass / mean_cyano_mass)
                total_ecw_biomass = float(args.od.split(",")[1]) * Biomass2OD[1]
                n_ecw = round(total_ecw_biomass / mean_ecw_mass)
                InitialConditions["cyano"]["OD"] = float(args.od.split(",")[0])
                InitialConditions["ecw"]["OD"] = float(args.od.split(",")[1])
            else:
                n_cyanos = int(random.uniform(1, 100))
            n_ecw = 0
            n_cells = n_cyanos
            total_cyano_biomass = n_cyanos * mean_cyano_mass
            total_ecw_biomass = n_ecw * mean_ecw_mass
            InitialConditions["cyano"]["OD"] = total_cyano_biomass / Biomass2OD[0]
            InitialConditions["ecw"]["OD"] = total_ecw_biomass / Biomass2OD[1]
            cyGroup = "group CYANO type 1"
            ecwGroup = ""
            if args.division.lower() == "on":
                cyDiv = f"fix d1 CYANO divide {args.biodt} v_EPSdens v_divDia1 {random.randint(1,1e6)}"
            else:
                cyDiv = ""
            ecwDiv = ""
            masses = "c_myMass[1]"
        elif args.culture_type == "ax-e":
            cell_types = ["ecw"]
            if args.cells is not None:
                n_ecw = int(args.cells.split(",")[1])
            elif args.od is not None:
                total_cyano_biomass = float(args.od.split(",")[0]) * Biomass2OD[0]
                n_cyanos = round(total_cyano_biomass / mean_cyano_mass)
                total_ecw_biomass = float(args.od.split(",")[1]) * Biomass2OD[1]
                n_ecw = round(total_ecw_biomass / mean_ecw_mass)
                InitialConditions["cyano"]["OD"] = float(args.od.split(",")[0])
                InitialConditions["ecw"]["OD"] = float(args.od.split(",")[1])
            else:
                n_ecw = int(random.uniform(1, 100))
            n_cyanos = 0
            n_cells = n_ecw
            total_cyano_biomass = n_cyanos * mean_cyano_mass
            total_ecw_biomass = n_ecw * mean_ecw_mass
            InitialConditions["cyano"]["OD"] = total_cyano_biomass / Biomass2OD[0]
            InitialConditions["ecw"]["OD"] = total_ecw_biomass / Biomass2OD[1]
            cyGroup = ""
            ecwGroup = "group ECW type 1"
            cyDiv = ""
            if args.division.lower() == "on":
                ecwDiv = f"fix d2 ECW divide {args.biodt} v_EPSdens v_divDia2 {random.randint(1,1e6)}"
            else:
                ecwDiv = ""
            masses = "c_myMass[1]"
        InitialConditions["cyano"]["StartingCells"] = n_cyanos
        InitialConditions["ecw"]["StartingCells"] = n_ecw
        InitialConditions["cyano"]["initial_biomass"] = total_cyano_biomass
        InitialConditions["ecw"]["initial_biomass"] = total_ecw_biomass
        _logger.info(f"Cyanos = {n_cyanos}")
        _logger.info(f"Ecw = {n_ecw}")
        _logger.info(f"Cyano biomass = {total_cyano_biomass} kg")
        _logger.info(f"Ecw biomass = {total_ecw_biomass} kg")
        _logger.info(f'Starting Cyano OD {InitialConditions["cyano"]["OD"]}')
        _logger.info(f'Starting Ecw OD {InitialConditions["ecw"]["OD"]}')
        RUN_DIR = Path(
            f"runs/Run_{n_cyanos}_{n_ecw}_{IPTG:.2e}_{args.reps}_{today}_{random.randint(1,1e6)}"
        )
        if not os.path.isdir(RUN_DIR):
            os.mkdir(RUN_DIR)
        # TODO embed cell type into metadata file and generate cell type programmatically

        grids = int(args.grid)
        while True:
            if (
                InitialConditions["Dimensions"][0] * 1e6 % grids == 0
                and InitialConditions["Dimensions"][1] * 1e6 % grids == 0
                and InitialConditions["Dimensions"][2] * 1e6 % grids == 0
            ):
                Mesh = f'{int(InitialConditions["Dimensions"][0]*1e6/grids)} {int(InitialConditions["Dimensions"][1]*1e6/grids)} {int(InitialConditions["Dimensions"][2]*1e6/grids)}'
                break
            else:
                grids += 1

        NutesNum = len(InitialConditions["Nutrients"]["Concentration"])
        for r in range(1, int(args.reps) + 1):
            L = [
                " NUFEB Simulation\r\n\n",
                f"     {n_cells} atoms \n",
                f"     {len(cell_types)} atom types \n",
                f"     {NutesNum} nutrients \n\n",
                f'  0.0e-4   {InitialConditions["Dimensions"][0] :.2e}  xlo xhi \n',
                f'  0.0e-4   {InitialConditions["Dimensions"][1] :.2e}  ylo yhi \n',
                f'  0.0e-4   {InitialConditions["Dimensions"][2] :.2e}  zlo zhi \n\n',
                " Atoms \n\n",
            ]

            j = 1
            for c, CellType in enumerate(cell_types, start=1):
                sizes = np.random.uniform(
                    InitialConditions[CellType]["min_size"],
                    InitialConditions[CellType]["max_size"],
                    size=(InitialConditions[CellType]["StartingCells"],),
                )
                initial_biomass = InitialConditions[CellType]["initial_biomass"]
                if initial_biomass > 0:
                    while not (
                        0.999 * initial_biomass
                        < sum(sizes ** 3 / 6 * np.pi * 370)
                        < 1.001 * initial_biomass
                    ):
                        sizes = np.random.uniform(
                            InitialConditions[CellType]["min_size"],
                            InitialConditions[CellType]["max_size"],
                            size=(InitialConditions[CellType]["StartingCells"],),
                        )
                    n = 0
                    for i in range(j, InitialConditions[CellType]["StartingCells"] + j):
                        size = sizes[n]
                        x = random.uniform(
                            0 + size, InitialConditions["Dimensions"][0] - size
                        )
                        y = random.uniform(
                            0 + size, InitialConditions["Dimensions"][1] - size
                        )
                        z = random.uniform(
                            0 + size, InitialConditions["Dimensions"][2] - size
                        )
                        L.append(
                            f'     %d {c} {size :.2e}  {InitialConditions[CellType]["Density"]} {x :.2e} {y :.2e} {z :.2e} {size :.2e} \n'
                            % (i)
                        )
                        j += 1
                        n += 1

            L.append("\n")
            L.append(" Nutrients \n\n")
            for i, nute in enumerate(
                InitialConditions["Nutrients"]["Concentration"].keys()
            ):
                L.append(
                    f'     %d {nute} {InitialConditions["Nutrients"]["State"][nute]} {InitialConditions["Nutrients"]["xbc"][nute]} {InitialConditions["Nutrients"]["ybc"][nute]} {InitialConditions["Nutrients"]["zbc"][nute]} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} {InitialConditions["Nutrients"]["Concentration"][nute] :.2e} \n'
                    % (i + 1)
                )

            L.append("\n")
            L.append(" Type Name \n\n")
            for c, CellType in enumerate(cell_types, start=1):
                L.append(f"     {c} {CellType}  \n")
            L.append("\n")
            L.append(" Diffusion Coeffs \n\n")
            for key in InitialConditions["Diff_c"].keys():
                L.append(f'     {key} {InitialConditions["Diff_c"][key]} \n')
            L.append("\n")
            L.append(" Growth Rate \n\n")
            for CellType in cell_types:
                L.append(
                    f'     {CellType} {InitialConditions[CellType]["GrowthRate"]} \n'
                )
            L.append("\n")
            L.append(" Ks \n\n")
            for CellType in cell_types:
                k = f"     {CellType}"
                for key in InitialConditions[CellType]["K_s"].keys():
                    k = k + " " + str(InitialConditions[CellType]["K_s"][key])
                k = k + f" \n"
                L.append(k)
            L.append("\n")
            for key in InitialConditions["cyano"]["GrowthParams"].keys():
                L.append(" " + key + f" \n\n")
                for CellType in cell_types:
                    L.append(
                        f'     {CellType} {InitialConditions[CellType]["GrowthParams"][key]} \n'
                    )
                L.append("\n")

            L.append("\n\n")

            # write atom definition file
            f = open(RUN_DIR / "atom.in", "w+")
            f.writelines(L)

        # write initial conditions json file
        dumpfile = open(RUN_DIR / "metadata.json", "w")
        json.dump(InitialConditions, dumpfile, indent=6)
        dumpfile.close()
        # write Inputscript
        # open the file
        atom_file_path = RUN_DIR.resolve() / "atom.in"
        if args.lammps == True:
            lammps = "dump    id all custom 100 output.lammmps id type diameter x y z"
        else:
            lammps = ""
        if args.hdf5 == True:
            hdf5 = (
                "dump    traj all bio/hdf5 100 trajectory.h5 id type radius x y z con"
            )
        else:
            hdf5 = ""
        if args.vtk == True:
            vtk = "dump    du1 all vtk 100 atom_*.vtu id type diameter x y z"
            grid = "dump    du2 all grid 100 grid_%_*.vti con"
            vtk_tarball = "true"
        else:
            vtk = ""
            grid = ""
            vtk_tarball = "false"
        if args.sucrose_halt is not None:
            suc_halt = f"fix h2 all halt {args.sucrose_halt} v_suc <= 1e-19"
        else:
            suc_halt = ""

        filein = open(TEMPLATES_DIR / "inputscript.txt")
        # read it
        src = Template(filein.read())
        # do the substitution
        result = src.safe_substitute(
            {
                "n": n,
                "SucRatio": SucRatio,
                "SucPct": SucPct,
                "n_cyanos": n_cyanos,
                "n_ecw": n_ecw,
                "Replicates": args.reps,
                "IPTG": f"{IPTG:.0e}",
                "Timesteps": args.timesteps,
                "date": today,
                "CYANOGroup": cyGroup,
                "ECWGroup": ecwGroup,
                "Zheight": InitialConditions["Dimensions"][2],
                "CYANODiv": cyDiv,
                "ECWDiv": ecwDiv,
                "GridMesh": f'{int(InitialConditions["Dimensions"][0]*1e6/int(args.grid))} {int(InitialConditions["Dimensions"][1]*1e6/int(args.grid))} {int(InitialConditions["Dimensions"][2]*1e6/int(args.grid))}',
                "lammps": lammps,
                "hdf5": hdf5,
                "vtk": vtk,
                "grid": grid,
                "masses": masses,
                "mass_max": f"{max_mass:.2e}",
                "DiffusionSteps": args.niter,
                "atom_file_path": atom_file_path,
                "biodt": args.biodt,
                "sucrose_halt": suc_halt,
            }
        )
        f = open(RUN_DIR / f"Inputscript.lammps", "w+")
        f.writelines(result)

        x = int(InitialConditions["Dimensions"][0] * 1e6)
        y = int(InitialConditions["Dimensions"][1] * 1e6)
        z = int(InitialConditions["Dimensions"][2] * 1e6)

        # write slurm script
        # open the file
        filein = open(TEMPLATES_DIR / "slurm.txt")

        # read it
        src = Template(filein.read())
        # do the substitution
        result = src.safe_substitute(
            {"job": f"NUFEB_{n}", "USER": args.user, "VTK": vtk_tarball}
        )
        f = open(f"NUFEB_{today}.slurm", "w+")
        f.writelines(result)
        # write local run script
        # open the file
        filein = open(TEMPLATES_DIR / "local.txt")

        # read it
        src = Template(filein.read())
        # do the substitution
        result = src.safe_substitute(
            {
                "n": n,
                "SucRatio": SucRatio,
                "SucPct": SucPct,
                "n_cyanos": n_cyanos,
                "n_ecw": n_ecw,
                "Reps": args.reps,
            }
        )
        f = open(RUN_DIR / f"local_{n_cyanos}_{n_ecw}_{SucPct}.sh", "w+")
        f.writelines(result)
        for arg in vars(args):
            _logger.info(f"{arg}={getattr(args, arg)}")

        _logger.info("Script ends here")


def run():
    """Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    """
    main(sys.argv[1:])


if __name__ == "__main__":

    run()
