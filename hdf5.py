"save rays in a hdf5 object"

import h5py
import numpy as np

def save(rayMatrix, path, name):
     with h5py.File(path+name+'.h5','w') as raySaving:
          raySaving.create_dataset('rays', data=rayMatrix)
          
def read(path, name):
     with h5py.File(path+name+'.h5','r') as raySaving:
          rayMatrix=raySaving.get('rays')
          rayMatrix=np.array(rayMatrix)

     return(rayMatrix)
     
