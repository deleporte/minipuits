exec open('base-functions.py')
exec open('integration.py')
exec open('level.py')
exec open('matrix.py')
exec open('plot-sphere.py')

def showGS(N):
    M = matrix(N)
    values, vectors = scipy.linalg.eigh(M)
    v = vectors.transpose()[0]
    plot(level(v,N))
    
def showmimic(N):
    M = mimic_matrix(N)
    values, vectors = scipy.linalg.eigh(M)
    v = vectors.transpose()[0]
    plot(level(v,N))
    
