import numpy as np
import numpy.random as nprd
import scipy.sparse as spa
import scipy.linalg

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt

import threading
import copy


def z(theta, phi):
    mod = 1./np.tan(np.pi/2.*theta)
    return mod*np.cos(np.pi*phi)+mod*np.sin(np.pi*phi)*1j

def plotsphere(function):
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

def func(data,N,z,w):
    return abs(sum(np.sqrt(scipy.special.binom(N,i%(N+1))*scipy.special.binom(N,i/(N+1)))*c*z**(i%(N+1))*w**(i/(N+1))/((1+abs(z)**2)*(1+abs(w)**2))**(0.5*N) for (i,c) in enumerate(data)))

def sumrow(function,k):
    Th = np.arange(0.01,1,.3)
    Ph = np.arange(0,2.01,.3)
    r = np.zeros(Ph.size)
    fichier=file('row'+str(k),'w')
    for l in range(Ph.size):
        for i in range(Th.size):
            for j in range(Th.size):
                r[l]+=function(z(Th[i],Ph[k]),z(Th[j],Ph[l]))*np.sin(Th[i]*np.pi)*np.sin(Th[j]*np.pi)/Th.size**2
    np.save(fichier,r)
    print "Row ", k+1, " finished"
    fichier.close()
    return 1

class sumThread(threading.Thread): #we want to return a value
    def __init__(self, function, k):
        threading.Thread.__init__(self)
        self.function = function
        self.k = k
        self._return = None
    def run(self):
        self._return = sumrow(self.function,self.k)
    def join(self):
        threading.Thread.join(self)
        return self._return
    
def process(function):
    Ph = np.arange(0,2.01,.3)
    threads=[]
    G = np.zeros((Ph.size,Ph.size))
    for k in range(Ph.size):
        threads.append(sumThread(function, k))
        threads[k].start()
    return 1
    
def digits(base,number):
    digits=[]
    N=int(number)
    while N > 0:
        digits.append(N-base*int(N/base))
        N /= base
    while len(digits)<6:
        digits.append(0)
    return digits
    
def showmu():
    M=np.diag([1.,1.,0.5,0.5]*4)
    L=[(0,6),(0,7),(1,2),(2,1),(2,3),(3,2),(4,5),(4,6),(5,4),(5,6),(6,4),(6,5),(6,7),(6,0),(7,6),(7,0)]
    aL=[(0, 1), (0, 2), (4, 2), (4, 3)]

    for l in L:
        M[2*l[0],2*l[1]]=0.25
        M[2*l[0]+1,2*l[1]+1]=-0.5

    mus=[]


    for t in range(10000):
        theta=np.pi*2*t/10000
        for l in aL:
            M[2*l[0],2*l[1]]=0.25*np.cos(theta)
            M[2*l[0]+1,2*l[1]+1]=-0.5*np.cos(theta)
            M[2*l[0],2*l[1]+1]=0.5*np.sin(theta)
            M[2*l[0]+1,2*l[1]]=0.25*np.sin(theta)
            M[2*l[1],2*l[0]]=M[2*l[0],2*l[1]]
            M[2*l[1]+1,2*l[0]+1]=M[2*l[0]+1,2*l[1]+1]
            M[2*l[1],2*l[0]+1]=M[2*l[0]+1,2*l[1]]
            M[2*l[1]+1,2*l[0]]=M[2*l[0],2*l[1]+1]
        vals=np.imag(scipy.linalg.eigvals(np.dot(J,M)))
        mu=0.
        for val in vals:
            if val>0:
               mu += val
        mus.append(mu)

def spin_matrices(N=6):
    Sz=spa.diags([np.array(range(N+1),dtype=float)-0.5*N],[0])
    vals=[]
    for k in range(N):
        m=-0.5*N+k
        vals.append(0.5*np.sqrt(0.5*N*(0.5*N+1)-m*(m+1)))
    Sx=spa.diags([vals,vals],[-1,1])
    Sy=spa.diags([np.array(vals)*1j,np.array(vals)*-1j],[-1,1])
    return [Sx,Sy,Sz]

def shift(N,mat,k):
    return spa.kron(spa.eye((N+1)**k),spa.kron(mat,spa.eye((N+1)**(4-k))))

class diagThread(threading.Thread): #we want to return a value
    def __init__(self, mat, vec):
        threading.Thread.__init__(self)
        self.mat = mat
        self.vec = vec
        self._return = None
    def run(self):
        self._return = iterate(self.mat,self.vec)
    def join(self):
        threading.Thread.join(self)
        return self._return

def iterate(mat,vec):
    for i in range(100):
        vec=mat.dot(vec)
        h=np.dot(np.conj(vec),vec)
        vec = vec/np.sqrt(h)
    return vec

def diagonalize(N=8):
    Sx=spin_matrices(N)[0]
    Sy=spin_matrices(N)[1]
    Sz=spin_matrices(N)[2]
    up=np.zeros((N+1))
    up[N]=1
    
    #downright=np.dot(scipy.linalg.expm(Sx.todense()*2.*np.pi*1j/3),up)
    upright=np.dot(scipy.linalg.expm(Sx.todense()*1.*np.pi*1j/3),up)
    upleft=np.dot(scipy.linalg.expm(-Sx.todense()*1.*np.pi*1j/3),up)
    
   #meanx=np.dot(np.conj(downright),np.array(np.dot(Sx.todense(),downright))[0])
   #meany=np.dot(np.conj(downright),np.array(np.dot(Sy.todense(),downright))[0])
   #meanz=np.dot(np.conj(downright),np.array(np.dot(Sz.todense(),downright))[0])
    meanxr=np.dot(np.conj(upright),np.array(np.dot(Sx.todense(),upright))[0])
    meanyr=np.dot(np.conj(upright),np.array(np.dot(Sy.todense(),upright))[0])
    meanzr=np.dot(np.conj(upright),np.array(np.dot(Sz.todense(),upright))[0])
    meanxl=np.dot(np.conj(upleft),np.array(np.dot(Sx.todense(),upleft))[0])
    meanyl=np.dot(np.conj(upleft),np.array(np.dot(Sy.todense(),upleft))[0])
    meanzl=np.dot(np.conj(upleft),np.array(np.dot(Sz.todense(),upleft))[0])
    print meanxr, meanxl
    print meanyr, meanyl
    print meanzr, meanzl
    temp=Sz*(meanzr)+Sx*meanxr+Sy*meanyr
    #Ham=shift(N,Sz*0.5*N,0)
    #Ham=Ham+shift(N,Sz*0.5*N,2)
    
    Ham=shift(N,temp,1)
    Ham=Ham+shift(N,temp,4)
    temp=Sz*(meanzl)+Sx*meanxl+Sy*meanyl
    Ham=Ham+shift(N,temp,0)
    Ham=Ham+shift(N,temp,2)
    
    Ham=Ham+shift(N,Sz,0).dot(shift(N,Sz,2))+shift(N,Sy,0).dot(shift(N,Sy,2))+shift(N,Sx,0).dot(shift(N,Sx,2))
    Ham=Ham+shift(N,Sz,0).dot(shift(N,Sz,1))+shift(N,Sy,0).dot(shift(N,Sy,1))+shift(N,Sx,0).dot(shift(N,Sx,1))
    Ham=Ham+shift(N,Sz,0).dot(shift(N,Sz,3))+shift(N,Sy,0).dot(shift(N,Sy,3))+shift(N,Sx,0).dot(shift(N,Sx,3))
    Ham=Ham+shift(N,Sz,1).dot(shift(N,Sz,3))+shift(N,Sy,1).dot(shift(N,Sy,3))+shift(N,Sx,1).dot(shift(N,Sx,3))
    Ham=Ham+shift(N,Sz,1).dot(shift(N,Sz,4))+shift(N,Sy,1).dot(shift(N,Sy,4))+shift(N,Sx,1).dot(shift(N,Sx,4))
    ##Finding the eigenvector with best eigenvalue using two threads

    Ham = Ham - 4*N*N*spa.eye((N+1)**5)
    #max_itr = (100.*np.log(N+1))**2
    v0=nprd.rand((N+1)**5*2)
    v0.dtype=complex
    v1=nprd.rand((N+1)**5*2)
    v1.dtype=complex
    cont=True
    k=0
    while cont:
        thread0 = diagThread(Ham,v0)
        thread1 = diagThread(Ham,v1)
        thread0.start()
        thread1.start()
        while thread0.isAlive() or thread1.isAlive():
            pass
        v0=thread0.join()
        v1=thread1.join()
        h01 = np.dot(np.conj(v0),v1)
        k+=100
        print k
        if 1.-abs(h01)<1e-5:
            cont=False
            print k
    return v0

    ## while cont:
    ##     thread.start_new_thread(iterate,(Ham,v0))
    ##     thread.start_new_thread(iterate,(Ham,v1))
    ## while cont:
    ##     v0 = Ham.dot(v0)
    ##     v1 = Ham.dot(v1)
    ##     h0 = np.dot(np.conj(v0),v0)
    ##     h1 = np.dot(np.conj(v1),v1)
    ##     v0 = v0/np.sqrt(h0)
    ##     v1 = v1/np.sqrt(h1)
    ##     h01 = np.dot(np.conj(v0),v1)
    ##     k+=1
    ##     #if k>max_itr:
    ##         #print "Warning: slow convergence"
    ##         #cont=False
    ##     if 1.-abs(h01)<1e-5:
    ##         cont=False
    ##         print k
    ## return v0
        
    #OP=Ham-sigma*spa.eye((N+1)**6)
    #return spa.linalg.LinearOperator(matvec=lambda v:spa.linalg.minres(OP,v,tol=1e-5)[0],shape=Ham.shape, dtype=Ham.dtype)
    #print "Bob"
   # return spa.linalg.eigsh(Ham,k=1,sigma=-20*N*N)[1].T[0]

def showmargins(N=8,n=3, vect=[0],j=0):
    mat=np.zeros((N+1,(N+1)**4),dtype=complex)
    for k in range((N+1)**5):
        lin=digits(N+1,k)[n]
        col=sum((N+1)**i*c for (i,c) in enumerate(digits(N+1,k)[:n]+digits(N+1,k)[n+1:]))
        mat[lin,col]=vect[k]
    (U,S,V)=scipy.linalg.svd(mat,full_matrices=False,overwrite_a=True)
    print "Coeff: ", S[j]**2
    Sx=spin_matrices(N)[0]
    #temp=np.dot(scipy.linalg.expm(Sx.todense()*2.*np.pi*1j/3),U.T[j])
    plotsphere(lambda z:abs(sum(np.sqrt(scipy.special.binom(N,i))*c*z**i/(1+abs(z)**2)**(0.5*N)for (i,c) in enumerate(U.T[j]))))

def findSVD(N=8,vect=[0]):
    mat=np.zeros(((N+1)**2,(N+1)**3),dtype=complex)
    for k in range((N+1)**5):
        lin=digits(N+1,k)[3]+(N+1)*digits(N+1,k)[4]
        col=sum((N+1)**i*c for (i,c) in enumerate(digits(N+1,k)[:3]))
        mat[lin,col]=vect[k]
    (U,S,V)=scipy.linalg.svd(mat,full_matrices=False,overwrite_a=True)
    return (U,S,V)

def builddoublemargins(N=8,U=np.array([0]),j=0):
    data=U.T[j]
    data.shape=(N+1,N+1)
    Sx=spin_matrices(N)[0].todense()
    data=np.dot(data,scipy.linalg.expm(Sx*1.*np.pi*1j/3))
    data=np.dot(scipy.linalg.expm(-Sx*1.*np.pi*1j/3),data)
    data.shape=((N+1)**2)
    process(lambda z,w: func(data,N,z,w))


def showdoublemargins(mat):
    plt.imshow(mat,cmap='hot',interpolation='nearest')
    plt.colorbar()
    plt.show()
