#Imports
import os
import subprocess
from pathlib import Path
from nufeb_tools import utils,plot
import pandas as pd
from string import Template

def main():
    TEMPLATES_DIR = (Path(__file__).parent) / 'templates'
    #Define inputs

    D = 370
    alpha =.2
    beta = 4
    delta=0.03

    #Change input params
    NUFEB_DIR = Path('~/NUFEB/src/USER-NUFEB/')

    filein = open( TEMPLATES_DIR / 'fix_bio_kinetics_monod.txt' )
            #read it
    src = Template( filein.read() )
            #do the substitution
    result = src.safe_substitute({'alpha' : alpha, 'beta' : beta, 'delta' : delta,
                                        
                                        })
    f= open(NUFEB_DIR / f"fix_bio_kinetics_monod.cpp","w+")
    f.writelines(result)
    #Compile NUFEB
    os.system("cd ~/NUFEB/")
    os.system("./install.sh --enable-hdf5 --enable-vtk")
    #Run simulation

    seed = subprocess.run(["nufeb-seed", "--cells", "100,0", "--d", "1e-4,1e-4,1e-4", "--grid", "20", "--t", "10000"])
    print("The exit code was: %d" % seed.returncode)

    #Extract output
    BASE_DIR = Path(f'runs/')
    folders = [path for path in BASE_DIR.iterdir() if path.is_dir()]
    data = [utils.get_data(directory=str(x)) for x in folders]
    Volume = 1e-4*1e-4*1e-4 #m^3
    CellNum2OD = Volume*1e6/0.3e-8
    SucroseMW = 342.3
    dfs = []
    for x in data:
        temp = pd.concat([x.ntypes.cyano/CellNum2OD,x.ntypes.step/60/60*x.timestep,x.avg_con.Sucrose.reset_index(drop=True)/SucroseMW*1e3],axis=1)
        temp.columns=['OD750','Hours','Sucrose']
        temp['SucroseExport'] = x.sucRatio/100
        dfs.append(temp)
    df = pd.concat(dfs)
    if not df.empty():
        print('success')
    #Compare output with experimental data

    #Optimize
if __name__ == "__main__":

    main()