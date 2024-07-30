#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
# coding: utf-8
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
from scipy.stats import loguniform
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

class Morphometries():	
    morphometries = {
         'talik' : np.array([[0,6500],
                            [0.25,6450],
                            [0.5,6400],
                            [0.75,6350],
                            [1,6300],
                            [1.25,6150],
                            [1.5,6000],
                            [1.75,3000],
                            [2,2500],
                            [2.25,2400],
                            [2.5,2300],
                            [2.75,1],
                            ]),
        'deep_box' : np.array([[0,6500],
                            [0.25,6450],
                            [0.5,6400],
                            [0.75,6350],
                            [1,6300],
                            [1.25,6250],
                            [1.5,6200],
                            [1.75,6150],
                            [2,6100],
                            [2.25,6050],
                            [2.5,3000],
                            [2.75,1]]),
        
        'deep_u' : np.array([[0,6500],
                            [0.25,6450],
                            [0.5,6400],
                            [0.75,6350],
                            [1,6300],
                            [1.25,6150],
                            [1.5,6000],
                            [1.75,5700],
                            [2,5400],
                            [2.25,4900],
                            [2.5,3000],
                            [2.75,1]
                            ]),
        
        'shallow_box' : np.array([[0,6500],
                                [0.25,6450],
                                [0.5,6400],
                                [0.75,6350],
                                [1,6300],
                                [1.25,6000],
                                [1.5,1]
                                ]),

        'shallow_u' : np.array([[0,6500],
                                [0.25,6200],
                                [0.5,5900],
                                [0.75,5500],
                                [1,5000],
                                [1.25,4000],
                                [1.5,1]
                                ])

       
    }
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
    def __init__(self,work_dir = "../LAKE_2024_JAMES",path_to_model_dir = "../LAKE_2024_JAMES",project_name = "TKL873"):
        #change this to your location on your VMs
        self.project_name = project_name
        self.path_to_model_dir = path_to_model_dir
        self.file_setup = os.path.join(self.path_to_model_dir,f"setup/{self.project_name}_setup.dat")
        self.file_driver = os.path.join(self.path_to_model_dir,f"setup/{self.project_name}_driver.dat")
        self.file_data = os.path.join(self.path_to_model_dir,f"data/{self.project_name}.dat")
        
        self.file_meteo = os.path.join(self.path_to_model_dir,f"meteo/{self.project_name}.dat")
        self.number = 0
        self.target_string = ""

        self.file_list = []
        self.target_values = []
        self.list_of_modifies_files = []
        #change this parameter depending on your username in google cloud
        self.setup_folder = ""
        self.work_dir = work_dir
        self.morphometries = Morphometries().morphometries
        self.completed_runs = []
        self.plots_output = "plots_output.pdf"
        

    def create_directories_auto(self, number_directories):
        '''
        Creates ``number_directories`` of directories and copies required files into them
        
        number_directories : int
        '''
        for i in range(number_directories):
            #change this to where you want to create all directories
            new_directory = os.path.join(self.work_dir,f"LAKE{i}")
            self.create_directory(new_directory)
            self.copy_required_files(new_directory)

    def copy_required_files(self, new_directory):
        #Copy necessary files to the new directory change this to your location on 
        #your computer if you are running from root directory then just leave it like that
        source_files = [
            "driver_file.dat",
            "results",
            "setup_file.dat",
            "setup",
            "launch",
            "crproj",
            "lake.out",
            "data",
            "meteo"
        ]

        for source in source_files:
            self.copy_to_directory(source, new_directory)

    def ig_f(self, dir, files):
        return [f for f in files if os.path.isfile(os.path.join(dir, f))]
        
    def copy_to_directory(self, source, destination_directory):
        source_path = os.path.join(self.path_to_model_dir, source) if self.path_to_model_dir else source
#         print(source_path)
        destination_path = os.path.join(destination_directory, source)
        if os.path.exists(destination_path):
#             print(f"Skipping copying {source} to {destination_path} as it already exists.")
            return

        if os.path.isfile(source_path):
            shutil.copy2(source_path, destination_path)
        elif os.path.isdir(source_path):
            shutil.copytree(source_path, destination_path, ignore = self.ig_f)
            
            
    def create_files(self):
        self.number = len(self.list_of_modifies_files)
        for i, (setup_file, driver_file) in enumerate(self.list_of_modifies_files):
            setup_filename = f"setup/{self.project_name}{i}_setup.dat"
            driver_filename = f"setup/{self.project_name}{i}_driver.dat"

            self.write_dictionary_to_file(setup_file, setup_filename)
            self.write_dictionary_to_file(driver_file, driver_filename)

        return self.file_list
    def _copy_files(self,args):
        source_file,destination_file = args
        with open(source_file, "r") as setup_file, open(destination_file, "w") as new_setup_file:
            shutil.copyfileobj(setup_file, new_setup_file)
    def create_files_single(self, source_directory, destination_directory,experiment_name):
        if not os.path.exists(destination_directory):
            os.makedirs(destination_directory)

#         print(f"Copying files from {source_directory} to {destination_directory}")

        for filename in os.listdir(source_directory):
            source_path = os.path.join(source_directory, filename)
            destination_path = os.path.join(destination_directory, filename)

            if source_path == destination_path:
#                 print(f"Skipping copying {source_path} to {destination_path} as it's the same file.")
                continue

#             print(f"Copying {source_path} to {destination_path}")

            # If the item is a file, copy it to the destination directory
            if os.path.isfile(source_path):
                shutil.copy2(source_path, destination_path)

            # If the item is a directory, recursively call create_files_single
            elif os.path.isdir(source_path):
                self.create_files_single(source_path, destination_directory,experiment_name)

    def create_data_file(self):
        try:
            if not os.path.exists("data"):
                os.makedirs("data")

            for i, (setup_file, driver_file) in enumerate(self.list_of_modifies_files):
                data_filename = os.path.join("data", f"YKD{i}.dat")
                self._copy_files((self.file_data, data_filename))

            print("Data files created successfully.")
        except Exception as e:
            raise ValueError(f"Error creating data files: {e}")

        try:
                if not os.path.exists("meteo"):
                    os.makedirs("meteo")
    
                for i, (setup_file, driver_file) in enumerate(self.list_of_modifies_files):
                    meteo_filename = os.path.join("meteo", f"TKL873{i}.dat")
                    self._copy_files((self.file_meteo, meteo_filename))
    
                print("Meteo files created successfully.")
        except Exception as e:
            raise ValueError(f"Error creating data files: {e}")
    def generate_and_write_files(self, num_files):
        """
        Generate and write tuples of driver and setup files to the "setup" folder.

        Args:
        num_files (int): The number of tuples to generate.
        """
        if not os.path.exists(os.path.join(self.work_dir, self.setup_folder)):
            os.makedirs(os.path.join(self.work_dir, self.setup_folder))


        for i in range(num_files):
            setup_folder = os.path.join(self.work_dir,f"LAKE{i}/setup")
            data_folder = os.path.join(self.work_dir,f"LAKE{i}/data")
            meteo_folder = os.path.join(self.work_dir,f"LAKE{i}/meteo")
            setup_filename = os.path.join(setup_folder, f"{self.project_name}{i}_setup.dat")
            driver_filename = os.path.join(setup_folder, f"{self.project_name}{i}_driver.dat")
            data_filename = os.path.join(data_folder,f"{self.project_name}{i}.dat")

            meteo_filename = os.path.join(meteo_folder,f"{self.project_name}{i}.dat")
   

            # Copy the content from the original setup and driver files to the new files
            with open(self.file_setup, "r") as src_setup_file, open(setup_filename, "w") as dst_setup_file:
                for line in src_setup_file:
                    dst_setup_file.write(line.replace("\t", " "))  # Replace tabs with spaces

            with open(self.file_driver, "r") as src_driver_file, open(driver_filename, "w") as dst_driver_file:
                for line in src_driver_file:
                    dst_driver_file.write(line.replace("\t", " "))  # Replace tabs with spaces
                    
            with open(self.file_data, "r") as src_driver_file, open(data_filename, "w") as dst_driver_file:
                for line in src_driver_file:
                    dst_driver_file.write(line.replace("\t", " "))  # Replace tabs with spaces

            with open(self.file_meteo, "r") as src_driver_file, open(meteo_filename, "w") as dst_driver_file:
                for line in src_driver_file:
                    dst_driver_file.write(line.replace("\t", " "))  # Replace tabs with spaces
                    
            # Append the tuple of filenames to the list_of_modifies_files
            self.list_of_modifies_files.append((setup_filename, driver_filename))


    def find_target(self, targets, new_values,number,is_mult = False,inflows = False):
        """
        Search for the specified target string in the setup and driver files and change its value.

        Args:
            target (str): The string to search for in the setup and driver files.
            new_values (list): A list of new values corresponding to each file in self.list_of_modifies_files.
        """
        # Make a copy of the new_values list to avoid modifying the original list
        
        if number != len(self.list_of_modifies_files) and not is_mult:
            raise ValueError(f"Number of new values({len(new_values)}) must be equal to the number of files in self.list_of_modifies_files({len(self.list_of_modifies_files)}).")
        if targets == ["ngrid_out"]:
            for (setup_file, driver_file), target_value in zip(self.list_of_modifies_files, new_values):
              

                print(f"Current target: {targets}, value: {target_value}")

                self.process_file_and_update_value(setup_file, targets[0], target_value, is_mult,inflows)
                self.process_file_and_update_value(driver_file, targets[0], target_value, is_mult,inflows)
            


        else:
            for (setup_file, driver_file), target_value in zip(self.list_of_modifies_files, new_values):

                print(target_value)
                for target,value in zip(targets,target_value):
                    print("Target values is",target_value)
                    if target == "dataname":
                        new_value = new_values
                        self.process_file_and_update_value(driver_file, target,target_value)
                        continue


                    self.process_file_and_update_value(setup_file, target,value,is_mult,inflows)

                    self.process_file_and_update_value(driver_file, target,value,is_mult,inflows)

    def process_file_and_update_value(self, file_path, target, new_value,is_mult = False,inflows=False):
        """
        Process the specified file, search for the target string, and update its value.

        Args:
            file_path (str): The path of the file to process.
            target (str): The string to search for in the file.
            new_value (str): The new value to replace the found target string.
        """
        with open(file_path, 'r') as file:
            lines = file.readlines()

        with open(file_path, 'w') as file:
            target_section = False
            for line in lines:
                line = line.strip()

                if line.startswith("#"):
                    # Skip lines starting with #
                    continue

                elif line.startswith("end"):
                    # Stop processing further lines after encountering "end"
                    file.write(line + "\n")
                    break
                elif (line.split()[0] == "fileinflow" or line.split()[0] == "fileoutflow") and new_value is None:
                    print(inflows)
                    print(file_path)
                    file.write(f"#{line}\n")
                    continue
                elif target in line and not is_mult:
                    file.write(f"{target} {new_value}\n")
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
                else:
                    # Copy the line as is
                    file.write(line + "\n")

    def read_file_path(self, filename):
        # Read the file and return the content as a string
        with open(filename, 'r') as file:
            content = file.read().strip()
        return content
    @staticmethod
    def create_directory(directory):
        os.makedirs(directory, exist_ok=True)
        print(f"Directory created: {directory}")
    def create_project_parallel(self, project_name,project_directory):
        
        if os.name == "posix":
            print("Operating system: Linux or macOS")
            OS = "linux" if os.uname().sysname == "Linux" else "OSX"
        else:
            print("Unknown operating system")
            return

        directories = [
            f"{project_directory}/results/{project_name}/everystep",
            f"{project_directory}/results/{project_name}/netcdf",
            f"{project_directory}/results/{project_name}/time_series",
            f"{project_directory}/results/{project_name}/hourly",
            f"{project_directory}/results/{project_name}/monthly",
            f"{project_directory}/results/{project_name}/daily",
        ]

        with multiprocessing.Pool() as pool:
            pool.map(self.create_directory, directories)
        
        
        # Modify driver file
        setup_folder = os.path.join(self.work_dir, project_directory)
        driver_file_path = os.path.join(project_directory, "driver_file.dat")
        print("Driver_file_path is",driver_file_path)
        if OS == "linux":
            sed_command = f"sed -i '2d' {driver_file_path} && sed -i \"\\$a setup/{project_name}_driver.dat\" {driver_file_path}"
        elif OS == "OSX":
            sed_command = f"sed -i '' '2d' {driver_file_path} && sed -i '' '$ a\\setup/{project_name}_driver.dat' {driver_file_path}"
        os.system(sed_command)

        # Modify setup file
        setup_file_path = os.path.join(project_directory, "setup_file.dat")
        print(setup_file_path)
        if OS == "linux":
            sed_command = f"sed -i '2d' {setup_file_path} && sed -i \"\\$a setup/{project_name}_setup.dat\" {setup_file_path}"
        elif OS == "OSX":
            sed_command = f"sed -i '' '2d' {setup_file_path} && sed -i '' '$ a\\setup/{project_name}_setup.dat' {setup_file_path}"
        os.system(sed_command)

        # Check if necessary files exist
        def check_file(file_path):
            if not os.path.isfile(file_path):
                print(f"Warning: The file {file_path} does not exist")

        check_file(f"./setup/{project_name}_setup.dat")
        check_file(f"./setup/{project_name}_driver.dat")
        check_file(f"./data/{project_name}.dat")
        check_file(f"./meteo/{project_name}.dat")

        print("Project for LAKE model created successfully.")

        
    #create multiple projects
    def create_multiple_projects(self,project_names,project_directories):
        """Uses multiprocessing pool to get mutliple project creation done in parallel and for each project calls 
        create_project_parallel function"""
        processes = []
        print("for plotting",project_directories)
        for project_name,project_directory in zip(project_names,project_directories):
            print(project_directory)
            process = multiprocessing.Process(target=self.create_project_parallel, args=(project_name,project_directory))
            process.daemon = False  # Set daemon to False to avoid the AssertionError
            processes.append(process)

        for process in processes:
            process.start()

        for process in processes:
            process.join()

        print("All projects created successfully.")

    def run_experiment(self, experiment_name,rundirectory,project_directory):
        """Intermidiate step in running the model,specifies the location and output of the result files"""
        if experiment_name:
            os.makedirs(f"results/{experiment_name}", exist_ok=True)

        self.run_model(rundirectory,project_directory)
        if experiment_name:
            basename = experiment_name.rsplit("_", 1)[0]
            source_directory = os.path.join(f"{project_directory}/results/{basename}")
#             Create a timestamp directories
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            destination_directory_with_timestamp = f"results/{experiment_name}_{timestamp}"

            # Copy files from source_directory to the new directory with a timestamp
            self.create_files_single(source_directory, destination_directory_with_timestamp, experiment_name)

    def generate_samples_for_SA(self, p_name, p_initial, perturbation, logparams, N, bounds = None, seed=''):
        params = []  # Dictionary of parameters
        if bounds:
            for name, p_bounds in zip(p_name, bounds):
                print(p_bounds)
                params.append(dict(name=name, bounds=p_bounds))
        else:        
            for name, init in zip(p_name, p_initial):
                p_bounds = [init - (init * perturbation), init + (init * perturbation)]
                params.append(dict(name=name, bounds=p_bounds, initial=init))

        # Set random seed if provided
        if seed != '':
            np.random.seed(int(seed))

        l = np.random.uniform(size=(N, len(params)))

        # Generate bounds, based on specification in params list
        lows = np.array([p['bounds'][0] for p in params])
        highs = np.array([p['bounds'][1] for p in params])

        # Figure out the spread, or difference between bounds
        spreads = highs - lows

        # Generate the sample matrix
        sm = l * spreads + lows

        # Apply loguniform for small param values only
        if len(logparams) > 0:
            inum = 0
            for ilog, p in zip(logparams, params):
                if ilog:
                    print(p['bounds'][0])
                    print(p['bounds'][1])
                    sm[:, inum] = loguniform.rvs(p['bounds'][0], p['bounds'][1], size=N)
                inum += 1

        return sm
    def auto_plot(self,directories,project_names):
        #plot ice thickness,co2,phosph,POCL
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

                for directory, project_name in zip(directories, project_names):
                    path_to_folder = os.path.join(directory, f"results/{project_name}/time_series")
                    path_to_file = os.path.join(path_to_folder, varfile[var])

                    if not os.path.exists(path_to_file):
                        print(f"File {path_to_file} does not exist. Skipping...")
                        continue

                    parsed_file = au.parse_univariate_file(path_to_file, var)

                    # Filtering out invalid data
                    filtered_data = parsed_file.loc[(parsed_file[var] != -999)]

                    sns.lineplot(data=filtered_data, x='Date', y=var, label=f'{project_name}')
                    

                plt.title(f'Time Series of {var.replace("_", " ").title()}')
                plt.xlabel('Date')
                plt.ylabel(var.replace("_", " ").title())
                plt.legend()

                # Save the current plot to the PDF file
                pdf.savefig()
                plt.close()
    
    def run_inflows_outflows(self,values,number):
        #values are going to be a dictionary where each key will be about
        return_list_inflows = []
        return_list_outflows = []
        for i,value in enumerate(values):
            directory = os.path.join(self.work_dir, f"LAKE{i}")
            path_to_meteo_file = os.path.join(directory,f"meteo/{self.project_name}{i}.dat")
            meteo_df = pd.read_csv(path_to_meteo_file)
            num_lines = len(meteo_df)
            inflows_df = pd.DataFrame(columns=['Date', 'width', 'U', 'temp', 'sal', 'Ux', 'Uy', 'DOC', 'POC', 'DIC', 'CH4', 'O2'])

            for key in value.keys():
              
         
                print("Number of lines is",num_lines)

                inflows_df['Date'] = pd.to_datetime(meteo_df[['Year', 'Month', 'Day']]).dt.strftime('%d%m%Y')
                inflows_df = inflows_df.fillna(-999)
#                 inflows_df[key] = value.get(key, [-999])[0]
                inflows_df[key] = [value.get(key, [-999][0])] * num_lines

    #                 inflows_df['U'] = values.get('U', [((4.3 * 1e-2) / 86400)]*len(inflows_df))

            inflows_path = os.path.join(directory, f'meteo/{self.project_name}{i}inflows.dat')
            inflows_df.to_csv(inflows_path, index=False, sep=' ', header=False)
            outflows_df = pd.DataFrame(columns=['Date', 'width', 'U'])
            outflows_df['Date'] = inflows_df['Date']
            outflows_df['width'] = inflows_df['width']
            outflows_df['U'] = inflows_df['U']
            outflows_path = os.path.join(directory,f'meteo/{self.project_name}{i}outflows.dat')
            outflows_df.to_csv(outflows_path, index=False, sep=' ', header=False)
            return_list_inflows.append([f"'{self.project_name}{i}inflows.dat'"])
        
            return_list_outflows.append([f"'{self.project_name}{i}outflows.dat'"])
        number_of_inflows = len(return_list_inflows)
        difference = number - number_of_inflows 
        tribheat = []
        other = []
        N_tribin = []
        for i in range(number_of_inflows):
            tribheat.append([2])
            other.append([1])
            N_tribin.append([1])
        self.find_target(["tribheat"],tribheat,number)
#         self.find_target(["N_tribin"],N_tribin,number)
        self.find_target(["N_triblev"],other,number)
        self.find_target(["iefflloc"],other,number)
        
        for i in range(difference):
            return_list_inflows.append([None])
            return_list_outflows.append([None])
        
        self.find_target(["fileinflow"],return_list_inflows,number,inflows = True)
        self.find_target(["fileoutflow"],return_list_outflows,number,inflows = True)

        
    def run_experiment_parallel(self, experiment_names,rundirectory,project_directories):
        """Using multiprocessing pool runs all of the processes in parallel, 
        which in turn uses self.run_model to run each process"""
        processes = []
        for experiment_name,project_directory in zip(experiment_names,project_directories):
            process = multiprocessing.Process(target=self.run_experiment, args=(experiment_name, rundirectory,project_directory))
            processes.append(process)

        for process in processes:
            process.start()

        for process in processes:
            process.join()

        print("All experiments completed successfully.")

    def run_model(self,rundirectory,project_directory):
        """Runs the model for each single model seperately """
        
        program_path = "./lake.out"
        print(project_directory)
        run_d = os.path.join(rundirectory,f"{project_directory}")
        try:
            completed_process = subprocess.run(
                program_path,
                shell=True,
                check=True,
                cwd=run_d,
                stdout=subprocess.PIPE,  # Capture standard output
                stderr=subprocess.PIPE,  # Capture standard error
                text=True               # Return output as text
            )
            self.completed_runs.append(run_d);
            print("Model run completed successfully.")
            print("Standard Output:")
            print(completed_process.stdout)
          

        except subprocess.CalledProcessError as e:
            print(f"Error occurred during model run. Exit code: {e.returncode}")
            print("Standard Output:")
            print(e.stdout)
            print("Standard Error:")
            print(e.stderr)
            raise e

    def clear(self):
        return self.list_of_modifies_files.clear()
    def change_morphometry(self,number):
        """Changes morphometry parameter in each driver file, 
        also changes required h10(max depth) parameter in driver file.
        and changes ngrid_out paramater in a setup file which is list of depths in every 0.2 m of specified morphometry
        Input parameter number just tells us muber of projects
        """

       
        
        
        result_list = []
        max_values = []
        
        result_max = []
        return_ngrid_list = []
        for name, data in self.morphometries.items():
            #morphometry
            list_df = []
            df = pd.DataFrame(data)
            
            list_df.append(df)
            result_list.append(list_df)
            #h10
            max_value = np.abs(np.max(data.transpose()[0]))
            n_depths = np.floor_divide(max_value,0.2)+1
            result_max.append([max_value])
                               
            n_grid_list = []
            #ngrid_out
            curr = 0.2
            
            for i in range(int(n_depths)+1):
                if max_value < curr:
                    break
                else:
                    n_grid_list.append(curr)
                    curr+=0.2
            return_ngrid_list.append(n_grid_list)


            
            
              
       
        
        self.find_target(["morphometry"],result_list,number,is_mult = True)
        self.find_target(["h10"],result_max,5)
        self.find_target(["ngrid_out"],return_ngrid_list,number,is_mult=True)
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
            path_to_file = os.path.join(path_to_folder, varfile[var])

            fig, axes = plt.subplots(len(varfile), 1, figsize=(10, 2*len(varfile)), sharex=True, sharey=True)
            
            for i, (var, file_name) in enumerate(varfile.items()):
                path_to_file = os.path.join(path_to_folder, file_name)
                if os.path.exists(path_to_file):
                    run_data = pd.read_csv(path_to_file, delim_whitespace=True, header=None, names=['Date', var])
                    run_data['Date'] = pd.to_datetime(run_data['Date'])
                    
                    axes[i].plot(run_data['Date'], run_data[var], label=var)
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
            
            plt.savefig(os.path.join(directory, f'{project_name}_timeseries_plots.jpg'), dpi=300)
            plt.close(fig)
        
    def create_and_run_sensitivity_analysis(self, number, p_name, target_values, rundirectory):
        """Combines all of the methods needed to run the multiprocessing model"""
        directories = []
        project_names = []
        for i in range(number):
            directories.append(os.path.join(self.work_dir, f"LAKE{i}"))

            project_names.append(f"{self.project_name}{i}")
        print(directories)
        print(project_names)
        self.clear()
#         self.create_data_file()
        values = [{
           
            'U': (4.3 * 1e-2) / 86400,
             "width": 32   # Example values# Example values
                }]
        self.create_directories_auto(number)
        self.generate_and_write_files(number)
        self.find_target(["dataname"], project_names, number)
#         self.change_morphometry(number)
        self.run_inflows_outflows(values,number)
#         list_of_modified_files = self.find_target(p_name, target_values, number)
        self.create_multiple_projects(project_names, directories)
        self.run_experiment_parallel(project_names, rundirectory, directories)
        self.auto_plot(directories,project_names)
    
    def clear_workdir(self):
        folder_path = self.work_dir
        # Check if the folder exists
        if os.path.exists(folder_path):
            for filename in os.listdir(folder_path):
                file_path = os.path.join(folder_path, filename)

            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)  # Remove files and symbolic links
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)  # Remove sub-directories
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")

        # If the folder doesn't exist, the code does nothing
        #ashutil.rmtree(self.work_dir)
    
#make it easier to change results folder


# In[9]:


sensitivity = SensitivityAnalysis()
project_names = []
directories = []

for i in range(5):
    directories.append(os.path.join("../LAKE_2024_JAMES", f"LAKE{i}"))

    project_names.append(f"TKL873{i}")
sensitivity.plot_depths(project_names,directories)


# In[7]:


sensitivity = SensitivityAnalysis()
p_name = ['khsO2', 'r0methprod',"VmaxCH4aeroboxid","khsCH4"]
p_initial = [2.1e+3,6.e+2,1.15e-7,3.75e+10]
perturbation = 0.75
logparams = np.zeros(len(p_initial))
logparams[1] = 1
N = 4
seed = ''
samples = sensitivity.generate_samples_for_SA(p_name, p_initial, perturbation, logparams, N, seed)
target_values = samples.tolist()
print(target_values)

number = len(target_values)
# number = sum(len(target_values[i]) for i in range(len(target_values)))
print(number)


# In[ ]:


# sensitivity.clear()
# sensitivity.create_directories_auto(number)
# sensitivity.generate_and_write_files(number)
# filename = ["setup","driver"]
# directory  = ["LAKE0","LAKE1","LAKE2","LAKE3","LAKE4"]
# sensitivity.create_data_file()
# sensitivity.find_target(["dataname"],["YKD0","YKD1","YKD2","YKD3","YKD4"],number,["driver","setup"])
# list_of_modified_files = sensitivity.find_target(p_name, target_values,number,filename)
# folder = ["YKD0","YKD1","YKD2","YKD3","YKD4"]
# # sensitivity.create_project_parallel(folder)
# sensitivity.create_multiple_projects(folder,directory)
# sensitivity.run_experiment_parallel(folder, rundirectory,directory)
rundirectory = os.path.abspath(os.getcwd())
# sensitivity.clear_workdir()
sensitivity. create_and_run_sensitivity_analysis(number,p_name,target_values,rundirectory)
# sensitivity.clear_workdir()
#all results are going to be saved with in root results folder





# In[ ]:


pip install openpyxl


# In[2]:


cd ..


# In[29]:


cd ../LAKE_2024_JAMES/LAKE0/TKL8730


# In[ ]:




