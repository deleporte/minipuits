def func(v,N,z):
    if len(v)<N+1:
        v.append(0)
        return func(v,N,z)
    else:
        somme=0
        for (k,vk) in enumerate(v):
            somme += vk*e(N,k,np.real(z), np.imag(z))
        return somme

def level(v,N):
    return lambda z:abs(func(v,N,z))**2.
