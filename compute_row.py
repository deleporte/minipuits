import sys
import numpy as np
import scipy.special

#argv[0] should be N
#argv[1] should be the name of the file containing the data
#argv[2] should be the number of rows
#argv[3] should be the row

def z(theta, phi):
    mod = 1./np.tan(np.pi/2.*theta)
    return mod*np.cos(np.pi*phi)+mod*np.sin(np.pi*phi)*1j

def func(data,N,z,w):
    return abs(sum(np.sqrt(scipy.special.binom(N,i%(N+1))*scipy.special.binom(N,i/(N+1)))*c*z**(i%(N+1))*w**(i/(N+1))/((1+abs(z)**2)*(1+abs(w)**2))**(0.5*N) for (i,c) in enumerate(data)))

#if len(sys.argv)>4:
#    continue

N=int(sys.argv[1])
data=np.load(sys.argv[2])
numrows=int(sys.argv[3])
row=int(sys.argv[4])
step=2.01/numrows

function = lambda z,w:func(data,N,z,w)

Th = np.arange(0.01,1,step)
Ph = np.arange(0,2.01,step)
r = np.zeros(Ph.size)
fichier=file(sys.argv[2]+'row'+str(row),'w')
for l in range(Ph.size):
    for i in range(Th.size):
        for j in range(Th.size):
            r[l]+=function(z(Th[i],Ph[row]),z(Th[j],Ph[l]))*np.sin(Th[i]*np.pi)*np.sin(Th[j]*np.pi)/Th.size**2
np.save(fichier,r)
print "Row ", row+1, " finished"
fichier.close()
