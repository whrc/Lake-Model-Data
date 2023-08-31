#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import  multi_proc_utils as mp

sa = mp.SensitivityAnalysis(work_dir='/home/amullen/Lake-Model-Data/multi_processing/outputs',
                            path_to_model_dir = "/home/amullen/LAKE")
#VmaxCH4aeroboxid 1.15E-7
#khsCH4 3.75E-2
#khsO2 2.1E-2
#r0methprod 6.E-8
#kc0 5.8E-6
#mubeta0 1.63E-4
p_name = ['kc0', 'mubeta0']
p_initial = [5.8E-6, 1.63E-4]
bounds = [[1e-7, 1e-1], [1e-5, 1e-1]]
perturbation = 1000
logparams = np.zeros(len(p_initial))
logparams[0] = 1
logparams[1] = 1
N = 10
seed = 2023
samples = sa.generate_samples_for_SA(p_name, p_initial, perturbation, logparams, N, bounds=bounds, seed=seed)
sample_values = samples.tolist()
print('Samples:',sample_values)
print('Sample size N:',N)

sa.file_setup = '/home/amullen/Lake-Model-Data/setup/YKD-burned_setup.dat'
sa.file_driver = '/home/amullen/Lake-Model-Data/setup/YKD-burned_driver.dat'
sa.file_data = '/home/amullen/Lake-Model-Data/data/prepped/YKD-burned.dat'
sa.file_meteo = '/home/amullen/Lake-Model-Data/data/prepped/YKD-burned.dat'

model_dir = '/home/amullen/LAKE'#os.path.abspath(os.getcwd())
sa.clear_workdir()
sa.create_and_run_sensitivity_analysis(N,p_name,sample_values,model_dir)
#all results are going to be saved with in root results folder
