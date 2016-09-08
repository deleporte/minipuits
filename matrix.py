import scipy.linalg

def matrix(N):
    M=np.empty([N+1,N+1], dtype=complex)
    for k in range(N+1):
        for l in range(N+1):
            M[k,l]=I(N,k,l)
    return M

def mimic_matrix(N):
    D = np.array(range(N+1)) - 0.5*N
    A = np.diag(D)
    A = np.dot(A,A)
    B = np.identity(N+1)*N*5.
    for k in range(N+1):
        if k>0:
            l = k-1
            B[k,l]=0.5*np.sqrt(0.25*N*(N+2)-(-0.5*N+k)*(-0.5*N+l))
        if k<N:
            l = k+1
            B[k,l]=0.5*np.sqrt(0.25*N*(N+2)-(-0.5*N+k)*(-0.5*N+l))
    return np.dot(A,B)+np.dot(B,A)
