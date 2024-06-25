#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import  multi_proc_utils as mp
import pandas as pd

#sa = mp.SensitivityAnalysis(work_dir='/home/amullen/Lake-Model-Data/model_output/YKD-SA/YKD-burned-july-start-vmaxch4', path_to_model_dir = "/home/amullen/LAKE")
#sa = mp.SensitivityAnalysis(work_dir='/home/amullen/Lake-Model-Data/model_output/YKD-SA/YKD-burned-july-start-r0methprod', path_to_model_dir = "/home/amullen/LAKE")
sa = mp.SensitivityAnalysis(work_dir='/home/amullen/Lake-Model-Data/model_output/YKD-SA/YKD-burned-july-start-CH4', path_to_model_dir = "/home/amullen/LAKE")
#VmaxCH4aeroboxid 1.15E-7
#khsCH4 3.75E-2
#khsO2 2.1E-2
#r0methprod 6.E-8
#kc0 5.8E-6
#mubeta0 1.63E-4

p_name = ['kc0', 'mubeta0', 'VmaxCH4aeroboxid', 'khsCH4', 'khsO2', 'r0methprod']
#p_initial = [5.8E-6, 1.63E-4, 1.15E-7, 3.75E-2, 2.1E-2, 6.0E-8]
p_initial = [5.8E-6, 1.63E-4, 1.11E-5, 9.5E-3, 2.1E-2, 6.0E-8]

#bounds = [[3.588e-7, 6.944e-7], [1e-7, 4e-6], [5.0E-6, 1.5E-5], [9.0E-3, 9.9E-3], [1.0E-2, 3.0E-2], [1e-5, 3e-05]] #best range for burned, good range but ch4 production too high
#bounds = [[4.50e-07, 4.55e-07], [1.8e-6, 2.2e-6], [6.4e-10, 1.0e-5], [9.0E-3, 9.9E-3], [1.0E-2, 3.0E-2], [4.5e-07, 4.6e-07]] #best range for burned, good range, but ch4 oxidation may be low

#bounds = [[3.35e-07, 3.45e-07], [2e-5, 2.3e-5], [9.5e-6, 1.2E-5], [6.0E-3, 8.2E-3], [1.3E-2, 4.1E-2], [6e-6, 1e-05]] #best range for unburned, good range but ch4 production too high
#bounds = [[3.35e-07, 3.45e-07], [1e-5, 2.3e-5], [1.0E-11, 5.0E-10], [6.0E-3, 8.2E-3], [1.0E-8, 5.0E-7], [1e-7, 5e-07]] #best range for unburned, good range, but ch4 oxidation may be low

#CH4 SA
p_name = ['r0methprod','VmaxCH4aeroboxid']
p_initial = [4.534203E-7, 1.15E-7]
bounds = [[2.499E-8 , 2.501e-8 ], [1.0E-9, 1.0E-5]] 

#CH4 SA r0methprod
#p_name = ['r0methprod']
#p_initial = [4.534203E-7]
#bounds = [[1.0E-8 , 2.0e-7 ]]

#CH4 SA VmaxCH4aeroboxid
#p_name = ['VmaxCH4aeroboxid']
#p_initial = [1.15E-7]
#bounds = [[1.0E-9, 1.0E-6]]


perturbation = 0.1
logparams = np.ones(len(p_initial))

N = 30
seed = 2023
samples = sa.generate_samples_for_SA(p_name, p_initial, perturbation, logparams, N, bounds=bounds, seed=seed)

sample_matrix_df = pd.DataFrame()

for i, p in enumerate(p_name):
    sample_matrix_df[p] = samples[:,i]

sample_matrix_df.to_csv(os.path.join(sa.work_dir,'sample_matrix.csv'))

sample_values = samples.tolist()
print('Samples:',sample_values)
print('Sample size N:',N)

sa.file_setup = '/home/amullen/LAKE/setup/YKD-burned-july-start_setup.dat'
sa.file_driver = '/home/amullen/LAKE/setup/YKD-burned-july-start_driver.dat'
sa.file_data = '/home/amullen/LAKE/data/YKD-burned-july-start.dat'
sa.file_meteo = '/home/amullen/LAKE/data/YKD-burned-july-start.dat'

model_dir = '/home/amullen/LAKE'#os.path.abspath(os.getcwd())
sa.clear_workdir()
sample_matrix_df.to_csv(os.path.join(sa.work_dir,'sample_matrix.csv'))
sa.create_and_run_sensitivity_analysis(N,p_name,sample_values,model_dir)
#all results are going to be saved with in root results folder
