import functools as ft

import shutil
import multiprocessing
import ast
import subprocess
import re
import datetime
import os
import shlex
from functools import partial
import numpy as np
from scipy.stats import loguniform, uniform
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.dates as mdates
import glob
import argparse
import yaml

class Auxilary_scripts():
    def parse_univariate_file(self,filepath, variable_name):
   
    
            file=pd.read_csv(filepath, delimiter=r"\s+", skiprows=7, index_col=None, header=None)
            file.columns=['year', 'month', 'day', 'hour', 'integration_time', 'depth', variable_name]
            file['Date'] = pd.to_datetime({'Year': file['year'], 'Month': file['month'], 'Day':file['day']})

            return file
    def parse_layer_file(self,filepath):
    
        layers=pd.read_csv(filepath, delimiter=r"\s+", skiprows=19, index_col=None, header=None)
        layers.columns=['year', 'month', 'day', 'hour', 'integration_time', 'water layer thickness, m', 
                               'W mixed layer thickness, m', 'E mixed layer thickness, m', 'S mixed layer thickness, m', 'N mixed layer thickness, m',
                               'W lower layer thickness, m', 'E lower layer thickness, m', 'S lower layer thickness, m', 'N lower layer thickness, m',
                               'ice layer thickness,   m', 'snow layer thickness,  m', 'bottom ice thickness,  m', 'reservoir volume,  m**3', 'volume deficit (accumulated),  m**3']
        layers['Date'] = pd.to_datetime({'Year': layers['year'], 'Month': layers['month'], 'Day':layers['day']})
        layers['mean_mixed_layer_thickness'] = layers[['W mixed layer thickness, m', 'E mixed layer thickness, m', 'S mixed layer thickness, m', 'N mixed layer thickness, m']].mean(axis=1)
        layers['mean_lower_layer_thickness'] = layers[['W lower layer thickness, m', 'E lower layer thickness, m', 'S lower layer thickness, m', 'N lower layer thickness, m']].mean(axis=1)
        layers = layers.drop(columns = ['W mixed layer thickness, m', 'E mixed layer thickness, m', 'S mixed layer thickness, m', 'N mixed layer thickness, m',
                                        'W lower layer thickness, m', 'E lower layer thickness, m', 'S lower layer thickness, m', 'N lower layer thickness, m'])
        layers['run_name'] = filepath.split('/')[-2]

        return layers

class log_wrapper:
    def __init__(self, cmdline, tag):
        self.cmdline = cmdline
        self.tag = tag

    def __enter__(self):
        print(f"Running: {self.cmdline}")

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type:
            print(f"Error occurred in {self.tag}: {exc_value}")
        else:
            print(f"{self.tag} completed successfully.")

class SensitivityAnalysis():
    '''
    This function lets LAKE model runs in parallel
    '''
    def __init__(self, conf):     

        self.conf=conf
        self.work_dir = self.conf['work_dir']
        
        self.project_name = self.conf['project_name']
        self.path_to_model_dir = self.conf['path_to_model_dir']

        #self.file_setup = os.path.join(self.path_to_model_dir,f"setup/{self.project_name}_setup.dat")
        #self.file_driver = os.path.join(self.path_to_model_dir,f"setup/{self.project_name}_driver.dat")
        #self.file_data = os.path.join(self.path_to_model_dir,f"data/{self.project_name}.dat")
        #self.file_meteo = os.path.join(self.path_to_model_dir,f"meteo/{self.project_name}.dat")

        self.n_runs = self.conf['n_runs']
        self.run_names = []
        self.run_dirs=[]

        for i in range(self.n_runs):
            run_name=f"{self.project_name}-{i}"
            self.run_names.append(run_name)
            self.run_dirs.append(os.path.join(self.work_dir, run_name))

        self.target_string = ""
        self.file_list = []
        self.target_values = []
        self.setup_driver_list = []
        
        #change this parameter depending on your username in google cloud
        self.setup_folder = ""

        #self.completed_runs = []
        self.plots_output = os.path.join(self.work_dir,'plots.pdf')

    def copy_required_files(self, new_directory, new_name):
        #Copy necessary files to the new directory change this to your location on 
        #your computer if you are running from root directory then just leave it like that
        os.mkdir(os.path.join(new_directory, 'results'))
        os.mkdir(os.path.join(new_directory, 'setup'))
        #os.mkdir(os.path.join(new_directory, 'data'))
        os.mkdir(os.path.join(new_directory, 'meteo'))
                        
        source_files = [
            os.path.join(self.path_to_model_dir, 'driver_file.dat'),
            #os.path.join(self.path_to_model_dir, 'data', f'{self.project_name}.dat'),
            os.path.join(self.path_to_model_dir, 'meteo', f'{self.project_name}.dat'),
            os.path.join(self.path_to_model_dir, 'setup_file.dat'),
            os.path.join(self.path_to_model_dir, 'setup', f'{self.project_name}_setup.dat'),
            os.path.join(self.path_to_model_dir, 'setup', f'{self.project_name}_driver.dat'),
            os.path.join(self.path_to_model_dir, 'crproj'),
            os.path.join(self.path_to_model_dir, 'launch'),
            os.path.join(self.path_to_model_dir, 'lake.out'),
            os.path.join(self.path_to_model_dir, 'tools'),
            os.path.join(self.path_to_model_dir, 'source'),
            os.path.join(self.path_to_model_dir, 'objfiles'),
            os.path.join(self.path_to_model_dir, 'results/debug'),
            os.path.join(self.path_to_model_dir, 'results/control_point')
        ]

        dests = [
            new_directory,
            #os.path.join(new_directory, 'data', f'{new_name}.dat'),
            os.path.join(new_directory, 'meteo', f'{new_name}.dat'),
            new_directory,
            os.path.join(new_directory, 'setup', f'{new_name}_setup.dat'),
            os.path.join(new_directory, 'setup', f'{new_name}_driver.dat'),
            new_directory,
            new_directory,
            new_directory,
            os.path.join(new_directory, 'tools'),
            os.path.join(new_directory, 'source'),
            os.path.join(new_directory, 'objfiles'),
            os.path.join(new_directory, 'results/debug'),
            os.path.join(new_directory, 'results/control_point')
        ]

        for source, dest in zip(source_files, dests):
            try:
                if os.path.isfile(source):
                    shutil.copy2(source, dest)
                elif os.path.isdir(source):
                    shutil.copytree(source, dest)
                else:
                    print(f'{source} does not exist')
                    break
            except:
                print(f'failed to copy source file {source} to destination {dest}')
                break
    
    def find_target(self, targets, new_values, number, is_mult = False, inflows = False):
        """
        Search for the specified target string in the setup and driver files and change its value.

        Args:
            targets (str): The string to search for in the setup and driver files.
            new_values (list): A list of new values corresponding to each file in self.setup_driver_list.
            is_mult (bool): new_value has multiple values corresponding to multiple lines in setup/driver file
        """
        # Make a copy of the new_values list to avoid modifying the original list
        
        if number != len(self.setup_driver_list) and not is_mult:
            raise ValueError(f"Number of new values({len(new_values)}) must be equal to the number of files in self.setup_driver_list({len(self.setup_driver_list)}).")
        
        if targets == ["ngrid_out"]:
            for (setup_file, driver_file), target_value in zip(self.setup_driver_list, new_values):
              
                print(f"Current target: {targets}, value: {target_value}")

                self.process_file_and_update_value(setup_file, targets[0], target_value, is_mult,inflows)
                self.process_file_and_update_value(driver_file, targets[0], target_value, is_mult,inflows)

        else:
            for (setup_file, driver_file), target_value in zip(self.setup_driver_list, new_values):

                for target,value in zip(targets,target_value):
                    #print("Target values is",value)
                    if target == "dataname":
                        self.process_file_and_update_value(driver_file, target,target_value)
                        continue

                    self.process_file_and_update_value(setup_file, target,value,is_mult,inflows)

                    #self.process_file_and_update_value(driver_file, target,value,is_mult,inflows)

    def process_file_and_update_value(self, file_path, target, new_value, is_mult = False, inflows=False):
        """
        Process the specified file, search for the target string, and update its value.

        Args:
            file_path (str): The path of the file to process.
            target (str): The string to search for in the file.
            new_value (str): The new value to replace the found target string.
            is_mult (bool): new_value has multiple values corresponding to multiple lines in setup/driver file
        """
        with open(file_path, 'r') as file:
            lines = file.readlines()

        with open(file_path, 'w') as file:
            target_section = False
            value_updated = False

            for line in lines:
                line = line.strip()

                if line.startswith("#"):
                    # Skip lines starting with #
                    continue

                elif line.startswith("end"):
                    # Stop processing further lines after encountering "end"
                    if value_updated==False:
                        file.write(f'{target} {new_value}\n{line}\n')
                    else:
                        file.write(line + "\n")
                    break

                elif (line.split()[0] == "fileinflow" or line.split()[0] == "fileoutflow") and new_value is None:
                    file.write(f"#{line}\n")
                    value_updated = True
                    continue
                elif target in line and not is_mult:
                    file.write(f"{target} {new_value}\n")
                    value_updated = True
                elif target in line and is_mult:
                    if isinstance(new_value, list):

                        file.write(f"{target} {len(new_value)} \n")
                        file.write("\n".join(map(str, new_value)) + "\n")
                    else:

                        file.write(f"{target} {len(new_value)} \n")
                        file.write(new_value.to_string(index=False, header=False) + "\n")
                    target_section = True
                elif re.match(r'^\d', line) and target_section:
                    continue
                elif not re.match(r'^\d', line) and target_section:
                    target_section = False
                    file.write(line + "\n")
                    value_updated = True
                else:
                    # Copy the line as is
                    file.write(line + "\n")

    def gen_samples(self, bounds, log_param, n_runs):

        if log_param:
            samples = loguniform.rvs(bounds[0], bounds[1], size=n_runs)
        else:
            samples = uniform.rvs(bounds[0], bounds[1]-bounds[0], size=n_runs)
        return samples

    def perturb_params(self, p_names, N, p_initial=None, perturbations=None, logparams=None, seed=None):
        # Set random seed if provided
        if seed is not None:
            np.random.seed(int(seed))

        samples = []
                        
        for name, init, perturbation, log_param in zip(p_names, p_initial, perturbations, logparams):
                
            if isinstance(perturbation, (list, tuple)):
 
                if len(perturbation) == 2:
                    samples.append(self.gen_samples(perturbation, log_param, N))

                elif int(len(perturbation)) == int(N):
                    samples.append(perturbation)
                    
                else:
                    print('if passing a list as perturbation, must be a range of length 2,')
                    print('or explicit values for each run of length n_runs')

            if isinstance(perturbation, float):

                p_bounds = [init - (init * perturbation), init + (init * perturbation)]
                samples.append(self.gen_samples(p_bounds, log_param, N))
        
        samples = np.array(samples)

        return samples.transpose()
    
    def mod_inflows_outflows_steady_state(self, path_to_meteo_file, setup_file, inflow_keys, inflow_values, 
                                          years_to_apply=None, months_to_apply=None):
        
        #values are going to be a dictionary where each key will be about
        inflows_path = path_to_meteo_file.split('.dat')[0]+'_inflows.dat'
        outflows_path = path_to_meteo_file.split('.dat')[0]+'_outflows.dat'
        meteo_df = pd.read_csv(path_to_meteo_file)
        num_lines = len(meteo_df)
        inflows_df = pd.DataFrame(columns=['Date', 'width', 'U', 'temp', 'sal', 'Ux', 'Uy', 'DOC', 'POC', 'DIC', 'CH4', 'O2'])
        inflows_df['Date'] = pd.to_datetime(meteo_df[['Year', 'Month', 'Day']]).dt.strftime('%Y%m%d')
        inflows_df = inflows_df.fillna(-999)
        
        if years_to_apply is None:
            years_to_apply = meteo_df['Year'].unique()
        
        if months_to_apply is None:
            months_to_apply = meteo_df['Month'].unique()
        #populate inflow columns based on supplied dictionary
        for key, value in zip(inflow_keys, inflow_values): 
            inflows_df.loc[(pd.to_datetime(inflows_df['Date']).dt.year.isin(years_to_apply)) & 
                           (pd.to_datetime(inflows_df['Date']).dt.month.isin(months_to_apply)), key] = value
            
            inflows_df.loc[~(pd.to_datetime(inflows_df['Date']).dt.year.isin(years_to_apply)) | 
                           ~(pd.to_datetime(inflows_df['Date']).dt.month.isin(months_to_apply)), key] = 0

        inflows_df.to_csv(inflows_path, index=False, sep=' ', header=False)
        outflows_df = pd.DataFrame(columns=['Date', 'width', 'U'])
        outflows_df['Date'] = inflows_df['Date']
        outflows_df['width'] = inflows_df['width']
        outflows_df['U'] = inflows_df['U']
        outflows_df.to_csv(outflows_path, index=False, sep=' ', header=False)

        self.process_file_and_update_value(setup_file, target='tribheat', new_value=2)
        #self.process_file_and_update_value(setup_file, target='tribheat', new_value=0)
        #self.find_target(["N_tribin"],1,1)
        self.process_file_and_update_value(setup_file, target='N_triblev', new_value=1)
        self.process_file_and_update_value(setup_file, target='iefflloc', new_value=1)
        self.process_file_and_update_value(setup_file, target='fileinflow', new_value=inflows_path.split('/')[-1])
        self.process_file_and_update_value(setup_file, target='fileoutflow', new_value=outflows_path.split('/')[-1])
    
    def run_lake(self, run_name, run_dir):
        """Runs the model for each single model seperately """
        try:
            completed_process = subprocess.run(
                f"bash crproj {run_name}; ./lake.out",
                shell=True,
                check=True,
                cwd=run_dir,
                stdout=subprocess.PIPE,  # Capture standard output
                stderr=subprocess.PIPE,  # Capture standard error
                text=True               # Return output as text
            )
            print("Project created.")
            print("Standard Output:")
            print(completed_process.stdout)
          

        except subprocess.CalledProcessError as e:
            print(f"Error occurred during model run. Exit code: {e.returncode}")
            print("Standard Output:")
            print(e.stdout)
            print("Standard Error:")
            print(e.stderr)
            raise e
        
    def run_experiment_parallel(self):
        """Using multiprocessing pool runs all of the processes in parallel, 
        which in turn uses self.run_model to run each process"""

        processes = []
        for run_name, run_dir in zip(self.run_names, self.run_dirs):
            process = multiprocessing.Process(target=self.run_lake, args=(run_name, run_dir))
            processes.append(process)

        for process in processes:
            process.start()

        for process in processes:
            process.join()

        print("All experiments ran successfully.")

    def change_morphometry(self, morphometry_name, setup_file, driver_file, ngrid_out_interval = 0.2):
        """Changes morphometry parameter in driver file, 
        also changes required h10(max depth) parameter in driver file.
        and changes ngrid_out paramater in a setup file which is list of depths in every 0.2 m of specified morphometry
        Input parameter number just tells us muber of projects
        """
        morphometry = np.array(self.conf['morphometries'][morphometry_name])
        #h10
        max_depth = np.abs(np.max(morphometry.transpose()[0]))
        max_area = np.abs(np.max(morphometry.transpose()[1]))
        ngrid_out_depths = np.floor_divide(max_depth,ngrid_out_interval)+1
                 
        n_grid_list = []
        #ngrid_out
        curr = 0.0
            
        for i in range(int(ngrid_out_depths)+1):
            if curr >= max_depth:
                break
            else:
                n_grid_list.append(curr)
                curr+=ngrid_out_interval
        
        self.process_file_and_update_value(driver_file, target='morphometry', new_value=pd.DataFrame(morphometry), is_mult = True)
        self.process_file_and_update_value(driver_file, target='h10', new_value=max_depth)
        self.process_file_and_update_value(driver_file, target='area_lake', new_value=max_area)
        self.process_file_and_update_value(setup_file, target='ngrid_out', new_value=pd.DataFrame(n_grid_list), is_mult = True)

    def plot_depths(self,project_names,directories):
        varfile = {
            "t_water":"water_temp  1  1f2.dat",
            "ch4":"methane_water  1  1f2.dat",
            "sal_water":"sal_water  1  1f2.dat",
            "co2":"co2_water  1  1f2.dat",
            "do":"oxygen_water  1  1f2.dat",
            "phosph":"phosph_water  1  1f2.dat",
            "chla":"chla  1  1f2.dat",
            "pocl":"POCL  1  1f2.dat",
            "pocd":"POCD  1  1f2.dat",
            "doc":"DOC  1  1f2.dat",
            "sod":"sod  1  1f2.dat",
            "prodox":"prodox  1  1f2.dat"}
        
        for project_name, directory in zip(project_names, directories):
            path_to_folder = os.path.join(directory, f"results/{project_name}/time_series")
            path_to_file = os.path.join(path_to_folder, varfile['sod'])

            fig, axes = plt.subplots(len(varfile), 1, figsize=(10, 2*len(varfile)), sharex=True, sharey=True)
            
            for i, (var, file_name) in enumerate(varfile.items()):
                path_to_file = os.path.join(path_to_folder, file_name)
                if os.path.exists(path_to_file):
                    run_data = pd.read_csv(path_to_file, delim_whitespace=True, header=None, names=['Date', var])
                    run_data['Date'] = pd.to_datetime(run_data['Date'])
                    
                    axes[i].plot(run_data['Date'], run_data['sod'], label=var)
                    axes[i].set_facecolor('#d6d6d6')
                    axes[i].set_ylabel(var, fontsize=13)
                    axes[i].yaxis.set_label_position("right")
                    axes[i].tick_params(axis='both', which='major', labelsize=13)
                    axes[i].legend(loc='upper right')
            
            locator = mdates.MonthLocator()
            fmt = mdates.DateFormatter('%b')
            X = plt.gca().xaxis
            X.set_major_locator(locator)
            X.set_major_formatter(fmt)
            
            plt.xticks(rotation=0)
            plt.xlim(pd.to_datetime('2022-10-15'), pd.to_datetime('2023-05-15'))
            fig.supylabel('Measurements', fontsize=13)
            fig.tight_layout()
            
            #plt.savefig(os.path.join(directory, f'{project_name}_timeseries_plots.jpg'), dpi=300)
            plt.close(fig)

    def gen_bgc_sample_matrix(self):

        return self.perturb_params(self.conf['bgc_targets'], 
                       self.n_runs, p_initial=self.conf['bgc_init_values'], 
                       perturbations=self.conf['bgc_perturb'], 
                       logparams=self.conf['bgc_logparams'], 
                       seed=1)
        
    def gen_in_outflows_sample_matrix(self):
        
        return self.perturb_params(self.conf['inflow_targets'], 
                       self.n_runs, p_initial=self.conf['inflow_init_values'], 
                       perturbations=self.conf['inflow_perturbations'], 
                       logparams=self.conf['inflow_logparams'], 
                       seed=1)
        
    def gen_morphometry_sample_matrix(self):
        
        conf_morphometries = self.conf['morphometries']
        sa_morphometry_names = []

        if self.n_runs != len(conf_morphometries.keys()):
            print(f'mismatch between number of runs ({self.n_runs}) and number of morphometries')
            print(f'supplied in conf file ({len(conf_morphometries.keys())}), looping through morphometries to fill runs')
        
        j=0
        for i in range(0,self.n_runs):
            sa_morphometry_names.append(list(conf_morphometries.keys())[j])
            j+=1
            if j == len(conf_morphometries.keys()):
                j=0
        
        print(f'setting morphometries, using {sa_morphometry_names} from config file')
        
        return sa_morphometry_names
                       
    def create_sensitivity_analysis(self):
        """Combines all of the methods needed to run the multiprocessing model"""
        sample_matrix = pd.DataFrame()
        sample_matrix['sample'] = np.arange(self.n_runs)

        if self.conf['modules']['bgc']:
            self.bgc_samle_matrix = self.gen_bgc_sample_matrix()
            sample_matrix[self.conf['bgc_targets']] = self.bgc_samle_matrix

        if self.conf['modules']['morphometry']:
            self.sa_morphometry_names = self.gen_morphometry_sample_matrix()
            sample_matrix['morphometry_name'] = self.sa_morphometry_names

        if self.conf['modules']['inflows_outflows']:
            self.inflows_samle_matrix = self.gen_in_outflows_sample_matrix()
            sample_matrix[self.conf['inflow_targets']] = self.inflows_samle_matrix
        
        sample_matrix.to_csv(os.path.join(self.work_dir, 'sample_matrix.csv'), index=False)

        #generate working directories
        for i, (run_name, run_dir) in enumerate(zip(self.run_names, self.run_dirs)):

            setup_file = os.path.join(run_dir, 'setup', f'{run_name}_setup.dat')
            driver_file = os.path.join(run_dir, 'setup', f'{run_name}_driver.dat')
            meteo_file = os.path.join(run_dir, 'meteo', f'{run_name}.dat')

            #remove if already exists
            if os.path.isdir(run_dir):
                shutil.rmtree(run_dir)

            #create directories
            os.mkdir(run_dir)
            self.copy_required_files(run_dir, run_name)
            self.setup_driver_list.append((setup_file,driver_file))
            
            #change dataname in driver file to path to input data
            self.process_file_and_update_value(os.path.join(run_dir, 'setup', f'{run_name}_driver.dat'), 
                                               target='dataname', new_value=f'\'{run_name}\'')
           
            #set morphometries if necessary
            if self.conf['modules']['morphometry']:
                self.change_morphometry(self.sa_morphometry_names[i], setup_file, driver_file, ngrid_out_interval = 0.2)

            if self.conf['modules']['inflows_outflows']:
                self.mod_inflows_outflows_steady_state(meteo_file, setup_file, 
                                                       self.conf['inflow_targets'], self.inflows_samle_matrix[i],
                                                       self.conf['inflow_years'], self.conf['inflow_months'])
            
            if not self.conf['modules']['inflows_outflows']:
                self.process_file_and_update_value(setup_file, target='tribheat', new_value=0)
        
        #modify biogeochemical params if necesarry       
        if self.conf['modules']['bgc']:
            self.find_target(self.conf['bgc_targets'], self.bgc_samle_matrix, 
                             self.conf['n_runs'], is_mult = False, inflows = False)

    def remove_file_or_dir(self, file):
        if os.path.isfile(file):
            os.remove(file)
        elif os.path.isdir(file):
            shutil.rmtree(file)
        else:
            print(f'file not found: {file}')
            return

    def consolidate_outputs(self):
        for d, n in zip(self.run_dirs, self.run_names):

            #remove everything besides results folder
            for i in os.listdir(d):
                if 'results' not in i:
                    self.remove_file_or_dir(os.path.join(d, i))

            #move all results from time_series to root of run directory
            for i in os.listdir(os.path.join(d,'results', n, 'time_series')):

                shutil.move(os.path.join(d,'results', n, 'time_series', i), os.path.join(d,i))

            self.remove_file_or_dir(os.path.join(d,'results'))
    
    def auto_plot(self):
        
        """Plots all the specified parameters,uses varfile dictionary 
        for dynamic acccess of file names in the result folder, 
        all of the plots are saved in the location specified in this variable, 
        self.plots_output which is in pdf format
        Also, inputs are directories which is list of locations of 
        each project and also list of project_names"""
        
        au = Auxilary_scripts()
        plotting_vars = ['t_water','sal_water', 'ch4', 'co2', 'do', 'phosph', 'chla', 'pocl', 'pocd', 'doc', 'sod', 'prodox']

        varfile = {
            "t_water":"water_temp  1  1f2.dat",
            "ch4":"methane_water  1  1f2.dat",
            "sal_water":"sal_water  1  1f2.dat",
            "co2":"co2_water  1  1f2.dat",
            "do":"oxygen_water  1  1f2.dat",
            "phosph":"phosph_water  1  1f2.dat",
            "chla":"chla  1  1f2.dat",
            "pocl":"POCL  1  1f2.dat",
            "pocd":"POCD  1  1f2.dat",
            "doc":"DOC  1  1f2.dat",
            "sod":"sod  1  1f2.dat",
            "prodox":"prodox  1  1f2.dat"
        }
        with PdfPages(self.plots_output) as pdf:
            for var in plotting_vars:
                plt.figure(figsize=(12, 6))

                for run_dir, run_name in zip(self.run_dirs, self.run_names):
            
                    path_to_file = os.path.join(run_dir, varfile[var])

                    if not os.path.exists(path_to_file):
                        print(f"File {path_to_file} does not exist. Skipping...")
                        continue

                    parsed_file = au.parse_univariate_file(path_to_file, var)

                    if self.conf['plot_depth'] is not None:
                        file_depths = parsed_file['depth'].unique()
                        min_diff_index = np.argmin(np.abs(file_depths-self.conf['plot_depth']))
                        closest_depth = file_depths[min_diff_index]
        
                        parsed_file = parsed_file.loc[parsed_file['depth']==closest_depth]
                        depth_label = str(closest_depth)

                    else:
                        depth_label = 'depth avg.'
                    
                    # Filtering out invalid data
                    filtered_data = parsed_file.loc[(parsed_file[var] != -999)]

                    sns.lineplot(data=filtered_data, x='Date', y=var, label=f'{run_name}')
                    

                plt.title(f'Time Series of {var.replace("_", " ").title()}, plot depth: {depth_label}')
                plt.xlabel('Date')
                plt.ylabel(var.replace("_", " ").title())
                plt.legend()

                # Save the current plot to the PDF file
                pdf.savefig()
                plt.close()

    def run_sensitivity_analysis(self):
        print('running experiment in parallel')
        self.run_experiment_parallel()
        print('experiment completed, consolidating outputs')
        self.consolidate_outputs()
        print('plotting')
        self.auto_plot()

def main(config_file):
    # Read YAML file
    with open(config_file, 'r') as stream:
        conf = yaml.safe_load(stream)

    sensitivity = SensitivityAnalysis(conf)
    sensitivity.create_sensitivity_analysis()
    sensitivity.run_sensitivity_analysis()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
                    prog='multi.py',
                    description='Runs sensitivity analysis for the LAKE model',
                    epilog='usage: python multi.py path/to/yaml/config.yaml')
    parser.add_argument('config', help='path to .yaml config file for LAKE SA')

    args = parser.parse_args()

    main(args.config)
