import tkinter as tk
from tkinter import filedialog
import numpy as np
import os
import json
import logging
# import winsound
import time

def select_file():
    root = tk.Tk()
    root.withdraw()
    return filedialog.askopenfilename()

def select_folder(default_path = None):
    root = tk.Tk()
    root.withdraw()
    root.call('wm', 'attributes', '.', '-topmost', True)
    if default_path:
        return filedialog.askdirectory(parent=root, initialdir=default_path)
    return filedialog.askdirectory(parent=root)

def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    return folder_path

def save_to_json(file_path, name, data):
    # Recursive conversion to handle numpy arrays of numpy arrays
    def convert_numpy_to_list(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert numpy array to list
        elif isinstance(obj, list):
            return [convert_numpy_to_list(item) for item in obj]
        return obj
    data = convert_numpy_to_list(data)
    # if os.path.exists(file_path):
    #     with open(file_path, "r") as f:
    #         datafile = json.load(f)
    # else:
    #     datafile = {}
    # if name is not None:
    #     datafile[name] = data
    # else: # if the data is already a dictionnaru ready to be saved
    #     datafile = data
    with open(file_path, "w") as f:
        json.dump(data, f, indent=4)

def save_to_txt(file_path, data):
    with open(file_path, "w") as f:
        f.write(data)

def load_from_json(file_path, name):
    if not os.path.exists(file_path):
        logging.error(f'{file_path} does not exist')
        return
    else:
        with open(file_path, "r") as f:
            datafile = json.load(f)
        return datafile[name]
      
def normalize(x, a=0, b=1):
    return (x - np.min(x)) / (np.max(x) - np.min(x)) * (b - a) + a

def visibility(x):
    return (np.max(x) - np.min(x)) / (np.max(x) + np.min(x))

def convert_power(w_power):
    return 10 * np.log10(w_power * 1000)

def beep(frequency=3000, duration=500):
    # winsound.Beep(frequency, duration)
    pass

def record_powermeter_data(powermeters, num_measurements, time_interval=1):
    time_arr = []
    measured_power_arr = [[] for _ in powermeters] 
    for i in range(num_measurements):
        print('%.2f ' % (i/num_measurements*100) + '%', end='\r')
        for powermeter in range(len(powermeters)):
            measured_power_arr[powermeter].append(powermeters[powermeter].read()) 
        time_arr.append(i)
        time.sleep(time_interval)
    print('100.00 %')
    return np.asarray(time_arr)*time_interval, [np.asarray(measured_power_arr[powermeter]) for powermeter in range(len(powermeters))]
