import numpy as np

def proj(u,v):
	return u*np.dot(v,u)/np.dot(u,u)

def GS(V):
	U=np.copy(V)
	for i in range(1,V.shape[0]):
		for j in range(0,i):
			U[i] -= proj(U[j],V[i])
		U[i] /= np.sqrt( (U[i]**2).sum(axis=0))
	return U
