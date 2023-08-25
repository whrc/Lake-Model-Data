#!/usr/bin/env python
# coding: utf-8

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
    def __init__(self,setup_folder = "setup"):
        #change this to your location on your VMs
        self.file_setup = "setup/YKD-burned_setup.dat"
        self.file_driver = "setup/YKD-burned_driver.dat"
        self.file_data = "data/YKD-burned.dat"
        self.number = 0
        self.target_string = ""
        self.file_list = []
        self.target_values = []
        self.dictionaries = []
        self.list_of_modifies_files = []
        self.setup_folder = setup_folder

    def create_directories_auto(self, number_directories):
        # Create the directories and copy necessary files to them
        for i in range(number_directories):
            #change this to where you want to create all directories
            new_directory = f"LAKE{i}"
            self.create_directory(new_directory)
            self.copy_required_files(new_directory)

    def copy_required_files(self, new_directory):
        # Copy necessary files to the new directory
        #change this to your location on your computer if you are running from root directory then 
        #just leave it like that
        source_files = [
            "driver_file.dat",
            "results",
            "setup_file.dat",
            "setup",
            "launch",
            "crproj",
            "lake.out",
            "data"
        ]

        for source in source_files:

            self.copy_to_directory(source, new_directory)

    def copy_to_directory(self, source, destination_directory):
        source_path = os.path.join(source)
#         print(source_path)
        destination_path = os.path.join(destination_directory, source)
        if os.path.exists(destination_path):
#             print(f"Skipping copying {source} to {destination_path} as it already exists.")
            return

        if os.path.isfile(source_path):
            shutil.copy2(source_path, destination_path)
        elif os.path.isdir(source_path):
            shutil.copytree(source_path, destination_path)
            
            
    def create_files(self):
        self.number = len(self.list_of_modifies_files)
        for i, (setup_file, driver_file) in enumerate(self.list_of_modifies_files):
            setup_filename = f"../setup/YKD{i}_setup.dat"
            driver_filename = f"../setup/YKD{i}_driver.dat"

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
    def generate_and_write_files(self, num_files):
        """
        Generate and write tuples of driver and setup files to the "setup" folder.

        Args:
        num_files (int): The number of tuples to generate.
        """
        if not os.path.exists(self.setup_folder):
            os.makedirs(self.setup_folder)

        for i in range(num_files):
            setup_folder = f"LAKE{i}/setup"
            setup_filename = os.path.join(setup_folder, f"YKD{i}_setup.dat")
            driver_filename = os.path.join(setup_folder, f"YKD{i}_driver.dat")

            # Copy the content from the original setup and driver files to the new files
            with open(self.file_setup, "r") as src_setup_file, open(setup_filename, "w") as dst_setup_file:
                for line in src_setup_file:
                    dst_setup_file.write(line.replace("\t", " "))  # Replace tabs with spaces

            with open(self.file_driver, "r") as src_driver_file, open(driver_filename, "w") as dst_driver_file:
                for line in src_driver_file:
                    dst_driver_file.write(line.replace("\t", " "))  # Replace tabs with spaces

            # Append the tuple of filenames to the list_of_modifies_files
            self.list_of_modifies_files.append((setup_filename, driver_filename))


    def find_target(self, targets, new_values,number):
        """
        Search for the specified target string in the setup and driver files and change its value.

        Args:
            target (str): The string to search for in the setup and driver files.
            new_values (list): A list of new values corresponding to each file in self.list_of_modifies_files.
        """
        # Make a copy of the new_values list to avoid modifying the original list
        
        if number != len(self.list_of_modifies_files):
            raise ValueError(f"Number of new values({len(new_values)}) must be equal to the number of files in self.list_of_modifies_files({len(self.list_of_modifies_files)}).")
        for (setup_file, driver_file), target_value in zip(self.list_of_modifies_files, new_values):
            for target,value in zip(targets,target_value):
                if target == "dataname":
                    new_value = new_values
                    self.process_file_and_update_value(driver_file, target,target_value)
                    continue

                self.process_file_and_update_value(setup_file, target,value)

                self.process_file_and_update_value(driver_file, target,value)

    def process_file_and_update_value(self, file_path, target, new_value):
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
            for line in lines:
                line = line.strip()

                if line.startswith("#"):
                    # Skip lines starting with #
                    file.write(line + "\n")
                elif line.startswith("end"):
                    # Stop processing further lines after encountering "end"
                    file.write(line + "\n")
                    break
                elif target in line:
                    file.write(f"{target} {new_value}\n")
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
    # Determine the operating system
        
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
        setup_folder = f"/home/kgurbanov/{project_directory}"
        driver_file_path = os.path.join(setup_folder, "driver_file.dat")
        if OS == "linux":
            sed_command = f"sed -i '2d' {driver_file_path} && sed -i \"\\$a setup/{project_name}_driver.dat\" {driver_file_path}"
        elif OS == "OSX":
            sed_command = f"sed -i '' '2d' {driver_file_path} && sed -i '' '$ a\\setup/{project_name}_driver.dat' {driver_file_path}"
        os.system(sed_command)

        # Modify setup file
        setup_file_path = os.path.join(setup_folder, "setup_file.dat")
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

        print("Project for LAKE model created successfully.")

        
    #create multiple projects
    def create_multiple_projects(self,project_names,project_directories):
        processes = []
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
        if experiment_name:
            os.makedirs(f"results/{experiment_name}", exist_ok=True)

        self.run_model(rundirectory,project_directory)
        if experiment_name:
            basename = experiment_name.rsplit("_", 1)[0]
            source_directory = f"{project_directory}/results/{basename}"

#             Create a timestamp directories
            timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
            destination_directory_with_timestamp = f"results/{experiment_name}_{timestamp}"

            # Copy files from source_directory to the new directory with a timestamp
            self.create_files_single(source_directory, destination_directory_with_timestamp, experiment_name)

    def generate_samples_for_SA(self, p_name, p_initial, perturbation, logparams, N, seed=''):
        params = []  # Dictionary of parameters

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
                    sm[:, inum] = loguniform.rvs(p['bounds'][0], p['bounds'][1], size=N)
                inum += 1

        return sm

    def run_experiment_parallel(self, experiment_names,rundirectory,project_directories):
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
    
    def create_and_run_sensitivity_analysis(self, number, p_name, target_values, rundirectory):
        directories = []
        project_names = []
        for i in range(number):
            directories.append(f"LAKE{i}")
            project_names.append(f"YKD{i}")
        print(directories)
        print(project_names)
        self.clear()
        self.create_directories_auto(number)
        self.generate_and_write_files(number)
        self.create_data_file()
        self.find_target(["dataname"], project_names, number)
        list_of_modified_files = self.find_target(p_name, target_values, number)
        self.create_multiple_projects(project_names, directories)
        self.run_experiment_parallel(project_names, rundirectory, directories)

