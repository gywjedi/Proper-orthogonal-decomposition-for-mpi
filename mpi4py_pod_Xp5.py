# -*- coding: utf-8 -*-
"""
Created on Sun Oct 29 19:17:50 2023

@author: ygn using pod to process Xp
"""

import pandas as pd
import numpy as np
import glob
import os
import math
from mpi4py import MPI
import scipy.io

def Xp_mol(base,input_dir,outdir,fnames):
    for fi in fnames:
        df = pd.read_csv(base + input_dir + fi,skiprows=9,sep="\s+",header=None)
        df.columns = ['id','type','x','y','z','xu','yu','zu']
        df = df.sort_values(by = ['id'])
        mol_beads = 1620;
        split_df = []
        num_mol = 100
        vec_X_arr_list = list()
        vec_Y_arr_list = list()
        vec_Z_arr_list = list()
        for i_mol in range(1,num_mol+1):
            _start_ = 1 + (i_mol-1)*mol_beads
            _end_ = mol_beads * i_mol
            #split_df.append(df.iloc[i_mol*mol_beads:(i_mol+1)*mol_beads])
            mol_df = df.query(f'id >= {_start_} and id <= {_end_}')
            xu = mol_df.xu.values
            yu = mol_df.yu.values
            zu = mol_df.zu.values
            
            mode_px = "./phi/phix" + f"{i_mol}" + ".mat"
            mode_py = "./phi/phiy" + f"{i_mol}" + ".mat"
            mode_pz = "./phi/phiz" + f"{i_mol}" + ".mat"
            loaded_matx = scipy.io.loadmat(mode_px)
            loaded_maty = scipy.io.loadmat(mode_py)
            loaded_matz = scipy.io.loadmat(mode_pz)
            x_phi = loaded_matx["phi"]
            y_phi = loaded_maty["phi"]
            z_phi = loaded_matz["phi"]
            
            vec_X_arr = []
            vec_Y_arr = []
            vec_Z_arr = []
            for p in range(mol_beads):
                xp = np.inner(xu,np.transpose(x_phi)[p])
                yp = np.inner(yu,np.transpose(y_phi)[p])
                zp = np.inner(zu,np.transpose(z_phi)[p])
                vec_X_arr.append(xp)
                vec_Y_arr.append(yp)
                vec_Z_arr.append(zp)
                
            vec_X_arr_list.append(vec_X_arr)
            vec_Y_arr_list.append(vec_Y_arr)
            vec_Z_arr_list.append(vec_Z_arr)
    
        rouseX_df = pd.DataFrame(data = np.transpose(vec_X_arr_list))
        rouseY_df = pd.DataFrame(data = np.transpose(vec_Y_arr_list))
        rouseZ_df = pd.DataFrame(data = np.transpose(vec_Z_arr_list))
    
        rouseX_df.to_csv(base + outdir + "x_" + fi,header=None,index=None,sep=",")
        rouseY_df.to_csv(base + outdir + "y_" + fi,header=None,index=None,sep=",")
        rouseZ_df.to_csv(base + outdir + "z_" + fi,header=None,index=None,sep=",")
   
if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    output_dir = "Xp5/"
    input_dir = r"_coord_/"
    base = './'
    try:
        os.mkdir(base + output_dir)
    except:
        print("folder has existed")
    fnames = [os.path.basename(x) for x in glob.glob(base + input_dir + 'dump.*.txt')]
    fnames.sort()
    total_numbers = len(fnames)
    number_per_node = total_numbers // size
    remainder = total_numbers % size
    start = rank * number_per_node + min(rank, remainder)
    end = start + number_per_node + (1 if rank < remainder else 0)
    local_results = Xp_mol(base,input_dir,output_dir,fnames[start:end])
    if rank == 0:
        print("single node")
