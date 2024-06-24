
from DOE_Reader import read_csv

import sys

from Solver import Solver

from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
 

def get_configurations(path):    
    return read_csv(path)



config_file_path = sys.argv[1]

configurations = get_configurations(config_file_path)

for cfg in configurations:
    print("EXPERIMENT " + cfg['#experiment'])
    my_solver = Solver()
    results = my_solver.solve(cfg)
    
    x = []
    y = []
    z = []
    
    for r in results:
        x.append(r[4])
        y.append(r[5])
        z.append(r[6])
        
    
    
  
    
    
    fig = plt.figure(figsize = (10, 7))
    ax = plt.axes(projection ="3d")
 
    
    ax.scatter3D(x, y, z, color = "green")
    plt.title("simple 3D scatter plot")
    
    
    id_name = '_OUTPUT/PLOTS/' + cfg['#problem'] + '/' + cfg['#instance_name'] + '-' + cfg['#algorithm'] + '-' + cfg['#problem'] + '-' + cfg['#configuration_name']
    plt.savefig(id_name)
    
   


