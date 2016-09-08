from scipy import integrate
import ctypes
import numpy as np
import os

pathlib = os.getcwd() + '/base-functions.so'
lib = ctypes.CDLL(pathlib)
mytype = (ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double, ctypes.c_double)
lib.realIntegrand.restype = ctypes.c_double
lib.realIntegrand.argtypes = mytype
lib.imagIntegrand.restype = ctypes.c_double
lib.imagIntegrand.argtypes = mytype

#def integrand(N,k,l,x,y):
#    return e(N,k,x,y)*h(x,y)*np.conj(e(N,l,x,y))/(1+x*x+y*y)**2

def I(N,k,l):
    return integrate.nquad(lambda x, y: lib.realIntegrand(N,k,l,x,y), [[-np.inf,np.inf],[-np.inf,np.inf]])[0] + integrate.nquad(lambda x, y: lib.imagIntegrand(N,k,l,x,y), [[-np.inf,np.inf],[-np.inf,np.inf]])[0]*1j

