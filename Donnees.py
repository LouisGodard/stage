"""In this module we define all the datas needed to compute the ray tracing"""

from Classes import *
import numpy as np

from sympy import *
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, V = symbols('f g V', cls=Function)
init_printing()

source=Source((0,0,2))  

rayon=Rayon(source.position,(20000,0,4600))


###sound speed analytic profil###

a=-0.1 #s^-1
V0=340 #m.s^-1 
V=V0+a*z #to work with sympy
#V=V0+a*log(z/0.1+1)
fV=lambda x,y,z:V0+a*z #for computational work
#fV=lambda x,y,z:V0+a*np.log(z/0.1+1)

#V=1500*(1.0+0.00737*(2*(z-1300)/1300-1+exp(-2*(z-1300)/1300))) 
#fV=lambda x,y,z:1500*(1.0+0.00737*(2*(z-1300)/1300-1+np.exp(-2*(z-1300)/1300)))
#https://oalib-acoustics.org/AcousticsToolbox/manual/node8.html#SECTION00340000000000000000


###ground###
# y=0

###ray tracing###

delta=0.02 #with the new formula of delta,it needs to be approximatly between
#0,001 et 1 to have a good precision and also a moderate computing time
distance=700
height=10

nbRay=21 #define to the number of rays shooted
angleMax=15/180*np.pi #shooting rays at regular interval between -angleMax and
#angleMax (in Rad !)
#if nbRay = 1 then one ray is traced at angleMax


###Saving###

plotting=True #to choose if you want to plot the rays
savePlot=False #to save the graph as a png with a tranparent background
hdf5=False #to save the raysMatrix as hdf5 file
text=False #to save the rays shooting under a text file

###ray plotting###

nbTotPoint=200
lineColor="r"
lineWidth=0.2

###ray hdf5 saving

path=''
name='rays'
