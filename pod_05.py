# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 19:28:56 2024

@author: ygn
"""

import numpy as np
import scipy
import glob
import os
import pandas as pd
import scipy.io
from mpi4py import MPI


def mpi_fnames_mol(fnames,mol_arr,polymerweight):
    for mol in mol_arr:
        list_of_df = []
        for fi in fnames:
            df = pd.read_csv(path + fi, sep="\s+",skiprows=9,header=None)
            df.columns = ['atid','atype','x','y','z','xu','yu','zu']
            df = df.sort_values(by = ['atid'])
            __start__ = 1 + (mol-1) * polymerweight
            __end__ = mol* polymerweight
            df_type1 = df.query(f'atid >= {__start__} and atid <= {__end__}')
            list_of_df.append(df_type1.zu.values)
        us = pd.DataFrame(data = list_of_df)

        ns, nx = us.shape
        dx1 = 1
        u_mean = np.zeros(ns)
        x = np.arange(1,nx+dx1,dx1)

        for i in range(ns):
            for j in range(nx):
                u_mean[i] = u_mean[i] + (us[j][i])/nx
        
        up = np.zeros((nx,ns))

        for i in range(ns):
            for j in range(nx):
                up[j][i] = us[j][i]-u_mean[i]
        
        C = np.zeros((nx,nx))
        for i in range(nx):
            for j in range(nx):
                if i < j:
                    pass
                else:
                    inner_p = np.zeros(ns)
                    for k in range(ns):
                        inner_p[k] = up[i][k]*up[j][k]
                
                    C[i][j]= np.trapz(inner_p,dx=dx1)/ns
                    if i == j:
                        pass
                    else:
                        C[j][i] = C[i][j]
                
        lam,v = np.linalg.eig(C)
        np.savetxt("lamz.txt",lam)
        idx = lam.argsort()[::-1]
        lam = lam[idx]
        v = v[:,idx]

        phi = np.zeros((nx,nx))
        for i in range(nx):
            for k in range(nx):
                phi[k][i] = v[k][i]###
        
        norm = np.zeros(nx)
        inner_p = np.zeros(nx)
        for p in range(nx):
            for j in range(nx):
                inner_p[j] = phi[j][p]*phi[j][p]
            norm[p]=np.sqrt(np.trapz(inner_p,dx=dx1))
            for j in range(nx):
                phi[j][p]=phi[j][p]/norm[p]

        str_emode="phiz" + f"{mol}" + ".mat"
        scipy.io.savemat(str_emode,{'phi':phi})  
        
if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    path = r"./_coord_/"
    fnames = [os.path.basename(x) for x in glob.glob(path + "dump.*.txt")]
    fnames.sort()
    total_numbers_mol = 100
    polymer_weight = 1620
    mol_per_node = total_numbers_mol // size
    remainder = total_numbers_mol % size
    start_f = rank * mol_per_node + min(rank, remainder) + 1
    end_f = start_f + mol_per_node + (1 if rank < remainder else 0)
    mol_list = [_ for _ in range(start_f, end_f)]
    mpi_fnames_mol(fnames,mol_list,polymer_weight)
