"""In this module we define all the datas needed to compute the ray tracing"""

from Classes import *
import numpy as np


from sympy import *
x, y, z = symbols('x y z')
V = symbols('V', cls=Function)

source=Source((0,0,2))  


###ground###
#y=0

###ray tracing###

delta=0.02 #with the new formula of delta,it needs to be approximatly between
#0,001 et 1 to have a good precision and also a moderate computing time
distance=330
height=8

###ray shooting###
nbRay=1 #define to the number of rays shooted
angleMax=-4.5/180*np.pi*38/40 #shooting rays at regular interval between -angleMax and
#angleMax (in Rad !)
#if nbRay = 1 then one ray is traced at angleMax


###sound speed analytic profil###

a=0.1 #s^-1
V0=340 #m.s^-1

V=V0+a*z #to work analytically with sympy
fV=lambda x,y,z:V0+a*z #for computational work

#V=V0+a*log(z/0.1+1)
#fV=lambda x,y,z:V0+a*np.log(z/0.1+1)

#V=1500*(1.0+0.00737*(2*(z-1300)/1300-1+exp(-2*(z-1300)/1300))) 
#fV=lambda x,y,z:1500*(1.0+0.00737*(2*(z-1300)/1300-1+np.exp(-2*(z-1300)/1300)))
#https://oalib-acoustics.org/AcousticsToolbox/manual/node8.html#SECTION00340000000000000000

analytic=0

###Sound speed discontinuous profil 2D###
zProfil=np.linspace(0,100*height,10001)
f=lambda z:V0+a*z
cProfil=f(zProfil)
#speedProfil=np.vstack((zProfil,cProfil))



###Saving###

plotting=True #to choose if you want to plot the rays
savePlot=True #to save the graph as a png with a tranparent background
hdf5Saving=True #to save the raysMatrix as hdf5 file
text=False #to save the rays shooting under a text file

path='' #don't forget the last /
name='rays'

###ray plotting###

nbPoint=200 #the number minimum of points per ray we want to trace the ray
lineColor="r"
lineWidth=0.2
fmt='png'



