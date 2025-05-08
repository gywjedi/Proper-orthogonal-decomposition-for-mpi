#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 10:51:06 2023

@author: ygn
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import re
import matplotlib as mpl
from mpi4py import MPI


def pdecay(pmodes,base):
    for p in pmodes:
        for count, _ in enumerate(base):
            fnames = [os.path.basename(x)[2:] for x in glob.glob(_ + 'Xp5/x_dump.*.txt')]
            x_rouse_arr = []
            y_rouse_arr = []
            z_rouse_arr = []
            fnames.sort()
            for fi in fnames:
                df_x = np.loadtxt(_ + 'Xp5/x_' + fi,delimiter=",",dtype=float,skiprows=int(p-1),max_rows=1,unpack=True)
                df_y = np.loadtxt(_ + 'Xp5/y_' + fi,delimiter=",",dtype=float,skiprows=int(p-1),max_rows=1,unpack=True)
                df_z = np.loadtxt(_ + 'Xp5/z_' + fi,delimiter=",",dtype=float,skiprows=int(p-1),max_rows=1,unpack=True)

                x_rouse_arr.append(df_x)
                y_rouse_arr.append(df_y)
                z_rouse_arr.append(df_z)
                
            xp_dataframe = pd.DataFrame(x_rouse_arr)
            yp_dataframe = pd.DataFrame(y_rouse_arr)
            zp_dataframe = pd.DataFrame(z_rouse_arr)
            
            T, nchains = xp_dataframe.shape
            
            rouse_modep = []        
            dt_range = np.arange(1,int(T/2))

            for dt in dt_range:
                xp0 = xp_dataframe.loc[range(T-dt)]
                xpt2 = np.sum(xp0.values*xp0.values)
                xpt = xp_dataframe.loc[range(dt,T)]
                xptxptdt = np.sum(xp0.values*xpt.values) 
                
                yp0 = yp_dataframe.loc[range(T-dt)]
                ypt2 = np.sum(yp0.values*yp0.values)
                ypt = yp_dataframe.loc[range(dt,T)]
                yptyptdt = np.sum(yp0.values*ypt.values) 

                zp0 = zp_dataframe.loc[range(T-dt)]
                zpt2 = np.sum(zp0.values*zp0.values)
                zpt = zp_dataframe.loc[range(dt,T)]
                zptzptdt = np.sum(zp0.values*zpt.values) 

                verhaltnis_dt = (xptxptdt + yptyptdt + zptzptdt)/(xpt2 + ypt2 + zpt2)
                rouse_modep.append(verhaltnis_dt)
        
        mydict = {"time":dt_range*10000*0.01,"norm_autocorr":rouse_modep}
        mode_p = pd.DataFrame(data = mydict)
        outfile_fname = f"./Pmode{p}.txt"
        mode_p.to_csv(outfile_fname,header=None,index=None,sep=" ")

if __name__ == '__main__':
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    base = ["./"]
    pmodes = [_ for _ in range(1,21)]
    total_numbers = len(pmodes)
    mode_per_node = total_numbers // size
    remainder = total_numbers % size
    start_mode = rank * mode_per_node + min(rank, remainder)
    end_mode = start_mode + mode_per_node + (1 if rank < remainder else 0)
    local_res = pdecay(pmodes[start_mode: end_mode], base)
