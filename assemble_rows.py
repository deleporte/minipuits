import sys
import numpy as np
import scipy.special

#argv[0] should be N
#argv[1] should be the name of the file containing the data
#argv[2] should be the number of rows

#if len(sys.argv)>3:
#    continue

N=int(sys.argv[1])
data=np.load(sys.argv[2])
numrows=int(sys.argv[3])
fichier=file(sys.argv[2]+"show",'w')
mat=[]
for row in range(numrows):
    vec=np.load(sys.argv[2]+"row"+str(row),'r')
    mat.append(list(vec))
matri=np.array(mat)
matri=matri/matri.max()
print matri
np.save(fichier,matri)
fichier.close()
