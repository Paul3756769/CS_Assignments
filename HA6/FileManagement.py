import os


def check_id(Firstline,Experiment_id,Measurement_id):
  firstline = Firstline.strip()
  # Check if the first line matches the criterion.
  id_values = firstline.split(";")
  experiment_id = 0
  measurement_id = 0
  if len(id_values) < 3:
    return False
  try:
    experiment_id = int(id_values[1])
    measurement_id = int(id_values[2])
  except ValueError:
    return False
  output = (experiment_id == Experiment_id) and (measurement_id == Measurement_id)
  
  return output

def collect_data(relative_path,experiment_id,measurement_id): 
  current_directory = os.getcwd()
  data_directory = os.path.join(current_directory,relative_path)

  dat_files = [f for f in os.listdir(data_directory) if f.endswith(".dat")]

  matched_files = []
  first_line = ""

  for dat_file in dat_files:
    dat_path = os.path.join(relative_path,dat_file)
    with open(dat_path, "r") as f:
      first_line = f.readline()

    if check_id(first_line,experiment_id,measurement_id):
      matched_files.append(dat_path)

  return matched_files

import numpy as np
import pandas as pd


def read_data_from_file(path,experiment_id):
    
    viable_experiments = [5,6] # Experiments with supported output format

    if experiment_id not in viable_experiments:
        raise ValueError("Experiment id is not recognized. Please choose from the following: " + str(viable_experiments))

    params = dict()
    with open(path, "r") as data_file:
        #skip first line
        data_file.readline()
        
        # read parameters
        
        # get parameter names as a list
        param_names = data_file.readline().replace('\n','').split(";")
        # get parameter values as a list
        param_values = data_file.readline().replace('\n','').split(";")

        # add to dictionary
        for i in range(len(param_names)):
            value = param_values[i]
            if value == "" or value == "\n":
                continue
            try:
                float_value = float(value)
                params[param_names[i]] = float_value
            except ValueError:
                print("Parameter " + param_names[i] + " is not recognized as a float. Value: " + value)
                #params[param_names[i]] = value
    
    
    # Read actual timeseries data
    
    # get  number of current line
    df = pd.read_csv(path, skiprows=7, sep=";", header=0)

    return df, params