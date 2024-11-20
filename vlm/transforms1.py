import numpy as np

def getDCM(csys1,csys2):
	dcm=np.zeros((3,3))
	for i in range(0,3):
		for j in range(0,3):
			dcm[i][j]=np.dot(csys1[i],csys2[j])
	return dcm

def transformVector(v,csys1,csys2):
	dcm=getDCM(csys1,csys2)
	return np.dot(dcm,v)

csys1=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]])
csys2=np.array([[ 1.,  0.,  0.], [ 0.,  -1.,  0.],  [ 0.,  0.,  -1.]])

v=np.array([0.70710678118654757,0.70710678118654757,0.70710678118654757])
v2=transformVector(v,csys2,csys1)
