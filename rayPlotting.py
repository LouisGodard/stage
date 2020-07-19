"This module manage the plotting of our rays"
import numpy as np
import matplotlib.pyplot as plt

def plot(raysMatrix,distance,height,lineColor='r',lineWidth=0.2,savePlot=True,plotting=False):
    for i in range(np.size(raysMatrix)//2):
        abscisses=raysMatrix[0,i]
        ordinates=raysMatrix[1,i]
        plt.plot(abscisses,ordinates,lineColor,linewidth=lineWidth)    
    sol=[0,distance]
    plt.plot(sol,[0,0],"--")
    plt.ylim(0,height)
    plt.xlim(0,distance)
    if savePlot:
        plt.savefig("result",transparent=True)
    if plotting:
        plt.show()
