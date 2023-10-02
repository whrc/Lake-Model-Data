#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import  multi_proc_utils as mp
import pandas as pd

sa = mp.SensitivityAnalysis(work_dir='/home/amullen/Lake-Model-Data/multi_processing/outputs',
                            path_to_model_dir = "/home/amullen/LAKE")
#VmaxCH4aeroboxid 1.15E-7
#khsCH4 3.75E-2
#khsO2 2.1E-2
#r0methprod 6.E-8
#kc0 5.8E-6
#mubeta0 1.63E-4

p_name = ['kc0', 'mubeta0', 'VmaxCH4aeroboxid', 'khsCH4', 'khsO2', 'r0methprod']
p_initial = [5.8E-6, 1.63E-4, 1.15E-7, 3.75E-2, 2.1E-2, 6.0E-8]
bounds = [[1e-8, 1e-4], [1e-6, 1e-2], [1e-9, 1e-5], [1e-4, 0.5], [1e-4, 0.5], [1e-10, 1e-6]]
#bounds = [[3e-8, 5e-7], [1e-4, 9e-3], [1e-6, 1e-5], [0.01, 0.03], [0.01, 0.09], [1e-5, 1e-4]]
perturbation = 0.1
logparams = np.zeros(len(p_initial))
logparams[0] = 1
logparams[1] = 1
logparams[2] = 1
logparams[3] = 1
logparams[4] = 1
logparams[5] = 1
N = 100
seed = 2023
samples = sa.generate_samples_for_SA(p_name, p_initial, perturbation, logparams, N, bounds=bounds, seed=seed)

sample_matrix_df = pd.DataFrame()

for i, p in enumerate(p_name):
    sample_matrix_df[p] = samples[:,i]

sample_matrix_df.to_csv(os.path.join(sa.work_dir,'sample_matrix.csv'))

sample_values = samples.tolist()
print('Samples:',sample_values)
print('Sample size N:',N)

sa.file_setup = '/home/amullen/LAKE/setup/YKD-burned-5y_setup.dat'
sa.file_driver = '/home/amullen/LAKE/setup/YKD-burned-5y_driver.dat'
sa.file_data = '/home/amullen/LAKE/data/YKD-burned-5y.dat'
sa.file_meteo = '/home/amullen/LAKE/data/YKD-burned-5y.dat'

model_dir = '/home/amullen/LAKE'#os.path.abspath(os.getcwd())
sa.clear_workdir()
sample_matrix_df.to_csv(os.path.join(sa.work_dir,'sample_matrix.csv'))
sa.create_and_run_sensitivity_analysis(N,p_name,sample_values,model_dir)
#all results are going to be saved with in root results folder
