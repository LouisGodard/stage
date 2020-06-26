"""In this module we define all the datas needed to compute the ray tracing"""

from Classes import *
import numpy as np

from sympy import *
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, V = symbols('f g V', cls=Function)
init_printing()

source=Source((0,0,2))  #créer une fonction plus tard

rayon=Rayon(source.position,(100,0,-4))


#définition du profil des vitesses et aussi penser aux données mesurées
a=0.1 #s^-1
V0=340 #m.s^-1 
V=V0+a*z #to work with sympy
#V=V0+a*log(z/0.1+1)
fV=lambda x,y,z:V0+a*z #for computational work
#fV=lambda x,y,z:V0+a*np.log(z/0.1+1)


#définition du sol

#a voir plus tard pour l'instant le sol est la droite y=0

#ray tracing

delta=0.2 #with the new formula of delta,it needs to be approximatly between
#0,05 et 1 to have a good precision and also a moderate computing time
#delta=0.0000000005
distance=300

#ray plotting
nbtot=200
 
