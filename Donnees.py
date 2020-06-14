"""In this module we define all the datas needed to compute the ray tracing"""

from Classes import *

from sympy import *
x, y, z, t = symbols('x y z t')
k, m, n = symbols('k m n', integer=True)
f, g, V = symbols('f g V', cls=Function)
init_printing()

source=Source((0,0,2))  #créer une fonction plus tard

rayon=Rayon(source.position,(100,0,8.5))


#définition du profil des vitesses et aussi penser aux données mesurées
a=0.1 #s^-1
V0=340 #m.s^-1 
V=V0+a*z #to work with sympy
fV=lambda x,y,z:V0+a*z #for computational work


#définition du sol

#a voir plus tard pour l'instant le sol est la droite y=0
