import numpy as np
import scipy.special
def e(N,k,x,y):
    return (x+y*1j)**k/(1.+x*x+y*y)**(N/2.)*np.sqrt(scipy.special.binom(N,k))

def h(x,y):
    return (x*x+y*y-1.)**2.*(1.+x*x+y*y+x)/(1.+x*x+y*y)**3.
