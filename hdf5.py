"save rays in a hdf5 object"

import h5py
import numpy as np

def save(raysMatrix, path, name):
     with h5py.File(path+name+'.h5','w') as raySaving:
          dt = h5py.vlen_dtype(np.dtype(np.float64))
          #dt = h5py.vlen_dtype(np.dtype([('absc', np.int64), ('ordi', np.float64)]))
          raySaving.create_dataset('raysData', data=raysMatrix, shape=(2,raysMatrix.size//2), dtype=dt)

def read(path, name):
     with h5py.File(path+name+'.h5','r') as raySaving:
          raysMatrix=raySaving.get('raysData')
          raysMatrix=np.array(raysMatrix)
     return(raysMatrix)
     
