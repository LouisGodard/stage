from RayShooting import *
import subprocess as sb
import datetime as dt
from Donnees import *
import rayPlotting
import hdf5

sb.run(["powercfg","/X","standby-timeout-ac","0"])
#to avoid standby of the computer when the program is executed and the computer
#is plugged

print("{} rays between {} and {} degrees".format(nbRay,-angleMax,angleMax))

start=dt.datetime.now()#start is the precise moment where we begin the rayshooting

#rayShooting()#execute rayShooting with datas in Donnees
raysMatrix=rayShooting()    

duration=dt.datetime.now()-start

print("ray shooting done")

print("{} rays between {} and {} rad".format(nbRay,-angleMax,angleMax))

print("Execution time: {}".format(duration))

sb.run(["powercfg","/X","standby-timeout-ac","15"])#restablish the
#standby parameter
#the number is the time (in minutes) you want your computer goes on stanby
#when he is plugged

if plotting or savePlot:
    rayPlotting.plot(raysMatrix,distance,height,lineColor,lineWidth,savePlot,plotting)

if hdf5Saving==True:
    hdf5.save(raysMatrix, path, name)

#sb.run("pause",shell=True)
