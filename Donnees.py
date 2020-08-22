"""In this module we extract from a text file all the datas needed to
compute the ray tracing"""

from Classes import *
import numpy as np

#in the case of the use of an analytic speed profile
from sympy import *
x, y, z = symbols('x y z')
V = symbols('V', cls=Function) 


###Extraction of the datas###

datas={"Zsource":2,"delta":0.02,"distance":1000,"height":100,"nbRay":11,"angleMax":45/180*np.pi,\
      "plotting":True,"savePlot":True,"hdf5saving":True,"path":'',"name":'rays',\
      "nbPoint":200,"lineColor":'r',"lineWidth":0.2,"fmt":'png',"analytic":1,"caracteristic":0}
#We stock datas in a dictionnary object because it is particulary appropriated
#you can modified the values in the datas dictionary if you want to run
#the program without a text file
V=340+0.1*z

textfilePath=input("""Enter the path of your file
Example : C:/Documents/data.txt
If your text file is in the same folder only the name is needed : Data.txt
You can use datas defined directly in the code by entering nothing.
Path:
""")

zProfil=np.array([])
cProfil=np.array([])
#the dimension of zProfil and cProfil are the same and together describe a
#speed profil
if textfilePath!='':#means the user want the program to use a data text file
    with open(textfilePath,'r') as data:

        #first we extract all datas except speed profile
        for ligne in data:
            if ligne[0]!='#'and ligne.rstrip('\n\r')!='' :#rstrip is used to remove the \n and this end of all line
                temp=ligne.rstrip('\n\r').split('=')
                try:
                    datas[temp[0]]=float(temp[1])
                except ValueError:#means the data is a string like the format
                    datas[temp[0]]=temp[1][1:-1]#[1:-1] used to remove first and
                    #last caracteres which normally are quotes
                if temp[0]=='analytic':
                    break
                
        #then we extract the speed profile
        begin=0
        for ligne in data:#in the case of a text file the loop continue where
            #the last loop ended
            if datas['analytic']:
                if ligne.rstrip('\n\r')=='%':
                    break
                elif ligne[0]!='#'and ligne.rstrip('\n\r')!='':
                    V=eval(ligne.rstrip('\n\r'))#V is now a function of sympy
                    #defined by the user
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
#A possible improvment is to make the program ables to treat more complex ground

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
name=datas['name'] #no extension needed in the name

###ray plotting###

nbPoint=int(datas['nbPoint']) #the number minimum of points per ray we want to trace the ray
lineColor=datas['lineColor']
lineWidth=datas['lineWidth']
fmt=datas['fmt']

###rays caracteristics###
caracteristic=datas['caracteristic']#to choose if you want to calculate the caracteristics of the rays

###Speed Profile###
analytic=int(datas['analytic'])

#V=V0+a*log(z/0.1+1)
#fV=lambda x,y,z:V0+a*np.log(z/0.1+1)

#V=1500*(1.0+0.00737*(2*(z-1300)/1300-1+exp(-2*(z-1300)/1300))) 
#fV=lambda x,y,z:1500*(1.0+0.00737*(2*(z-1300)/1300-1+np.exp(-2*(z-1300)/1300)))
#https://oalib-acoustics.org/AcousticsToolbox/manual/node8.html#SECTION00340000000000000000
