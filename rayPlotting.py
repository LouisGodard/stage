"This module manage the plotting of our rays"
import numpy as np
import matplotlib.pyplot as plt


def plot(raysMatrix,distance,height,fmt='png',nbPoint=200,lineColor='r',lineWidth=0.2,name="result",savePlot=True,plotting=False):
    
    x0=raysMatrix[0,3]
    abscisses=np.array([])
    ordinates=np.array([])
    step=distance/nbPoint
    for a,b,c,x1,reflexion in raysMatrix[1:]:
        if np.isnan(b):#means we will browse a new ray
            plt.plot(abscisses,ordinates,lineColor,linewidth=lineWidth)
            abscisses=np.array([])
            ordinates=np.array([])
            x0=x1
        else:
            absc=np.arange(x0,x1,step)
            Z=np.poly1d([a,b,c])
            ordi=Z(absc)
            abscisses=np.hstack((abscisses,absc))
            ordinates=np.hstack((ordinates,ordi))
            x0=x1
    plt.plot(abscisses,ordinates,lineColor,linewidth=lineWidth)    
    sol=[0,distance]
    plt.plot(sol,[0,0],"--")
    plt.ylim(0,height)
    plt.xlim(0,distance)
    if savePlot:
        plt.savefig(name+'.'+fmt,transparent=True, format=fmt)
    if plotting:
        plt.show()
