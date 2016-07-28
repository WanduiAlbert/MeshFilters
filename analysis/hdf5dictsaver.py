import h5py
import numpy as np

def save_dict_to_hdf5(group, mydict):
    for key in mydict:
        
        group.create_dataset(name=key)