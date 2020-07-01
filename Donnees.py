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

a=1 #s^-1
V0=340 #m.s^-1 
V=V0+a*z #to work with sympy
V=V0+a*log(z/0.1+1)
fV=lambda x,y,z:V0+a*z #for computational work
fV=lambda x,y,z:V0+a*np.log(z/0.1+1)

#V=1500*(1.0+0.00737*(2*(z-1300)/1300-1+exp(-2*(z-1300)/1300))) 
#fV=lambda x,y,z:1500*(1.0+0.00737*(2*(z-1300)/1300-1+np.exp(-2*(z-1300)/1300)))
#https://oalib-acoustics.org/AcousticsToolbox/manual/node8.html#SECTION00340000000000000000


###ground###
# y=0

###ray tracing###

delta=0.2 #with the new formula of delta,it needs to be approximatly between
#0,01 et 1 to have a good precision and also a moderate computing time
distance=325
height=6.5

nbRay=3 #define to the number of rays shooted
angleMax=4.5/180*np.pi #shooting rays at regular interval between -angleMax and
#angleMax (in Rad !)
#if nbRay = 1 then one ray is traced at angleMax

###ray plotting###

nbTotPoint=200
lineColor="r"
lineWidth=0.2
 
