"save rays in a hdf5 object"

import h5py
import numpy as np

def save(raysMatrix,raysIndex,raysData, path, name):
     with h5py.File(path+name+'.h5','w') as raySaving:
          raySaving.create_dataset('rays', data=raysMatrix)
          raySaving.create_dataset('index', data=raysIndex)
          raySaving.create_dataset('datas', data=raysData)
          
def read(path, name):
     with h5py.File(path+name+'.h5','r') as raySaving:
          raysMatrix=raySaving.get('rays')
          raysMatrix=np.array(raysMatrix)

          raysIndex=raySaving.get('index')
          raysIndex=np.array(raysIndex)

          raysData=raySaving.get('datas')
          raysData=np.array(raysData)

     return(raysMatrix,raysIndex,raysData)
     
