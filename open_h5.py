from __future__ import print_function

import h5py as h5


filename = "MPS_GS_8chi_Heisenberg3x3.h5"

def scan_hdf5(path, recursive=True, tab_step=2):
    def scan_node(g, tabs=0):
        print(' ' * tabs, g.name)
        for k, v in g.items():
            if isinstance(v, h5.Dataset):
                print(' ' * tabs + ' ' * tab_step + ' -', v.name)
            elif isinstance(v, h5.Group) and recursive:
                scan_node(v, tabs=tabs + tab_step)
    with h5.File(path, 'r') as f:
        scan_node(f)

        
def scan_hdf52(path, recursive=True, tab_step=2):
    def scan_node(g, tabs=0):
        elems = []
        for k, v in g.items():
            if isinstance(v, h5.Dataset):
                elems.append(v.name)
            elif isinstance(v, h5.Group) and recursive:
                elems.append((v.name, scan_node(v, tabs=tabs + tab_step)))
        return elems
    with h5.File(path, 'r') as f:
        return scan_node(f)
    
import h5py

f = h5.File(filename)
for key in f.keys():
    print(key) #Names of the root level object names in HDF5 file - can be groups or datasets.
    print(type(f[key])) #   # returns as a numpy array
group = f[key]

#Checkout what keys are inside that group.
for key in group.keys():
    print(key)

print("FOLDER INDENTED")
scan_hdf5(filename)
#print("PATH LIKE")
#print(scan_hdf52("MPS_GS_8chi_Heisenberg3x3.h5"))

#https://myhdf5.hdfgroup.org/