#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import  multi_proc_utils as mp

sa = mp.SensitivityAnalysis(work_dir='/home/ejafarov/LAKE/Lake-Model-Data/multi_processing/outputs',
                            path_to_model_dir = "/home/ejafarov/LAKE")
p_name = ['khsO2', 'r0methprod',"VmaxCH4aeroboxid","khsCH4"]
p_initial = [2.1e+3,6.e+2,1.15e-7,3.75e+10]
perturbation = 0.75
logparams = np.zeros(len(p_initial))
logparams[1] = 1
N = len(p_name)
seed = ''
samples = sa.generate_samples_for_SA(p_name, p_initial, perturbation, logparams, N, seed)
sample_values = samples.tolist()
print('Samples:',sample_values)
print('Sample size N:',N)

sa.file_setup = '/home/ejafarov/LAKE/Lake-Model-Data/setup/YKD-unburned_setup.dat'
sa.file_driver = '/home/ejafarov/LAKE/Lake-Model-Data/setup/YKD-unburned_driver.dat'
sa.file_data = '/home/ejafarov/LAKE/Lake-Model-Data/data/prepped/YKD-unburned.dat'

rundirectory = '/home/ejafarov/LAKE'#os.path.abspath(os.getcwd())
sa.clear_workdir()
sa.create_and_run_sensitivity_analysis(N,p_name,sample_values,rundirectory)
#all results are going to be saved with in root results folder
