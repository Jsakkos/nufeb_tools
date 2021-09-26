#Imports
import os
from random import uniform
import subprocess
from pathlib import Path
from nufeb_tools import utils,plot
import pandas as pd
from string import Template
import numpy as np
from hyperopt import fmin, tpe, hp,space_eval
from hyperopt.pyll import scope
from numpy.ma import MaskedArray
import sklearn.utils.fixes

sklearn.utils.fixes.MaskedArray = MaskedArray
from skopt import gp_minimize
import pickle
def main(x):

    alpha = x[0]
    beta = x[1]
    delta = x[2]
    mu = x[3]
    rho = x[4]
    exp_low = [1.38,.041872]
    exp_high = [1.146667,1.141355]
    TEMPLATES_DIR = (Path(__file__).parent) / 'templates'
    #Define inputs


    #Change input params
    HOME = Path.home()
    NUFEB_DIR = HOME / '/NUFEB/src/USER-NUFEB/'

    filein = open( TEMPLATES_DIR / 'fix_bio_kinetics_monod.txt' )
            #read it
    src = Template( filein.read() )
            #do the substitution
    result = src.safe_substitute({'alpha' : alpha, 'beta' : beta, 'delta' : delta,
                                        
                                        })
    os.chdir(HOME)
    f= open("NUFEB/src/USER-NUFEB/fix_bio_kinetics_monod.cpp","w+")
    f.writelines(result)
    #Compile NUFEB
    os.chdir(str(HOME / 'NUFEB'))
    os.system("./install.sh --enable-hdf5 --enable-vtk")
    #Clean old simulations
    os.system('nufeb-clean')
    #Run simulation
    text = f'nufeb-seed --cells 100,0 --d 1e-4,1e-4,1e-4 --grid 20 --t 10000 --mucya {mu} --sucR 0 --rhocya {rho}'
    os.system(text)
    text = f'nufeb-seed --cells 100,0 --d 1e-4,1e-4,1e-4 --grid 20 --t 10000 --mucya {mu} --sucR 1 --rhocya {rho}'
    os.system(text)
    os.system('./run.sh')
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
    low_suc = df.loc[(df.Hours > 23.8) & (df.Hours < 24) & (df.SucroseExport==0)].mean()[['OD750','Sucrose']].to_numpy()
    high_suc = df.loc[(df.Hours > 23.8) & (df.Hours < 24) & (df.SucroseExport==1)].mean()[['OD750','Sucrose']].to_numpy()

    #Compare output with experimental data
    return np.sqrt((low_suc - exp_low)**2).sum() + np.sqrt((high_suc - exp_high)**2).sum()
    #Optimize
def optimize():
    

    res = gp_minimize(main,                  # the function to minimize
                    [(0.1, .5),(1,5),(0.01,.1),(1e-5,1e-6),(320,390)],      # the bounds on each dimension of x
                    acq_func="EI",      # the acquisition function
                    n_calls=15,         # the number of evaluations of f
                    n_random_starts=5,  # the number of random initialization points
                    random_state=1234)
    print(res)
    file_pi = open('results.obj', 'w') 
    pickle.dump(res, file_pi)
    """         rho = 370
        alpha =.2
        beta = 4
        delta=0.03
        mu = 1.67e-5 """
"""     space = scope.main(hp.uniform('alpha', .01, .5),hp.uniform('beta', 1, 5),hp.uniform('delta', .01, .1),hp.uniform('mu', 1e-6, 1e-5),hp.uniform('rho',320,390))
    best = fmin(main,
    space=space,
    algo=tpe.suggest,
    max_evals=100)

    print(best)
    print(space_eval(space, best)) """
    
if __name__ == "__main__":

    optimize()