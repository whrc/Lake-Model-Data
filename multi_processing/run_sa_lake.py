#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import  multi_proc_utils as mpu 

p_name = ['khsO2', 'r0methprod',"VmaxCH4aeroboxid","khsCH4"]
p_initial = [2.1e+3,6.e+2,1.15e-7,3.75e+10]
logparams = [0,1,1,0]
variance = 0.75
N = len(p_initial)
seed = ''

sa = mpu.SensitivityAnalysis()
samples = sa.generate_samples_for_SA(p_name, p_initial, variance, logparams, N, seed)
sample_values = samples.tolist()
print('sample values:',sample_values)

number = len(sample_values)
print('number of samples:',number)

rundirectory = os.path.abspath(os.getcwd())
print('working dir:',rundirectory)
sa.create_and_run_sensitivity_analysis(number,p_name,sample_values,rundirectory)
