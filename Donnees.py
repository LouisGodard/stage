"""In this module we define all the datas needed to compute the ray tracing"""

from Classes import *
import numpy as np


from sympy import *
x, y, z = symbols('x y z')
V = symbols('V', cls=Function)


###Datas extraction from text file###

datas={"Zsource":2,"delta":0.02,"distance":1000,"height":100,"nbRay":11,"angleMax":45/180*np.pi,\
      "plotting":True,"savePlot":True,"hdf5saving":True,"path":'',"name":'rays',\
      "nbPoint":200,"lineColor":'r',"lineWidth":0.2,"fmt":'png'}
print("""Entrez le chemin d'accès à votre fichier
Exemple : C:/Documents/data.txt
Ne rentrez rien pour utiliser les données prédéfinies""")
textfilePath=input()
begin=0
zProfil=np.array([])
cProfil=np.array([])
if textfilePath!='':
    with open(textfilePath,'r') as data:
        for ligne in data:
            if ligne[0]!='#'and ligne.rstrip('\n\r')!='' :
                temp=ligne.rstrip('\n\r').split('=')
                try:
                    datas[temp[0]]=float(temp[1])
                except ValueError:
                    datas[temp[0]]=temp[1][1:-1]
                if temp[0]=='analytic':
                    break
        for ligne in data:
            if datas['analytic']:
                if ligne.rstrip('\n\r')=='%':
                    break
                elif ligne[0]!='#'and ligne.rstrip('\n\r')!='':
                    V=eval(ligne.rstrip('\n\r'))
            else:
                if begin==1 and ligne[0]!='#'and ligne.rstrip('\n\r')!='':
                    temp=ligne.rstrip('\n\r').split(':')
                    Ze=np.array([np.float(temp[0])])
                    Ve=np.array([np.float(temp[1])])
                    zProfil=np.concatenate((zProfil,Ze))
                    cProfil=np.concatenate((cProfil,Ve))
                if ligne.rstrip('\n\r')=='%':
                    begin=1
                    
print(datas)

source=Source((0,datas['Zsource']))  


###ground###
#y=0

###ray tracing###

delta=datas['delta'] #with the new formula of delta,it needs to be approximatly between
#0,001 et 1 to have a good precision and also a moderate computing time
distance=datas['distance']
height=datas['height']

###ray shooting###
nbRay=int(datas['nbRay']) #define to the number of rays shooted
angleMax=datas['angleMax']/180*np.pi #shooting rays at regular interval between -angleMax and
#angleMax (in Rad !)
#if nbRay = 1 then one ray is traced at angleMax



###Saving###

plotting=datas['plotting'] #to choose if you want to plot the rays
savePlot=datas['savePlot'] #to save the graph as a png with a tranparent background
hdf5Saving=datas['hdf5Saving'] #to save the raysMatrix as hdf5 file

path=datas['path'] #don't forget the last /
name=datas['name']

###ray plotting###

nbPoint=int(datas['nbPoint']) #the number minimum of points per ray we want to trace the ray
lineColor=datas['lineColor']
lineWidth=datas['lineWidth']
fmt=datas['fmt']


###sound speed analytic profil###

a=0.1 #s^-1
V0=340 #m.s^-1

#V=V0+a*z #to work analytically with sympy
if datas['analytic']:
    fV=lambdify((x,z),V,"numpy")


#V=V0+a*log(z/0.1+1)
#fV=lambda x,y,z:V0+a*np.log(z/0.1+1)

#V=1500*(1.0+0.00737*(2*(z-1300)/1300-1+exp(-2*(z-1300)/1300))) 
#fV=lambda x,y,z:1500*(1.0+0.00737*(2*(z-1300)/1300-1+np.exp(-2*(z-1300)/1300)))
#https://oalib-acoustics.org/AcousticsToolbox/manual/node8.html#SECTION00340000000000000000

analytic=int(datas['analytic'])

###Sound speed discontinuous profil 2D###
#zProfil=np.linspace(0,100*height,10001)
#f=lambda z:V0+a*z
#cProfil=f(zProfil)
#speedProfil=np.vstack((zProfil,cProfil))




