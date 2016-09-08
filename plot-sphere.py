from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

def z(theta, phi):
    mod = 1./np.tan(np.pi/2.*theta)
    return mod*np.cos(np.pi*phi)+mod*np.sin(np.pi*phi)*1j

def plot(function):
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.view_init(20,200)
    #X = np.arange(-0.8,0.8,.1)
    #Y = np.arange(-0.8,0.8,.1)
    #R = np.arange(0,1,.099)
    #T = np.arange(0,1.05,.05)
    #R, T = np.meshgrid(R,T)
    Th = np.arange(0.01,1,.01)
    Ph = np.arange(0,2.01,.01)
    Th, Ph = np.meshgrid(Th,Ph)
    Z = np.cos(np.pi*Th)
    X = np.sin(np.pi*Th)*np.cos(np.pi*Ph)
    Y = np.sin(np.pi*Th)*np.sin(np.pi*Ph)
    F = function(z(Th, Ph))
    N = F/F.max()
    surf = ax.plot_surface(X,Y,Z, rstride=1, cstride=1, facecolors=cm.jet(N), linewidth=0, antialiased=False, shade=False)
    plt.show()
