import numpy as np

def getNearestCPData(cpData,zPos):
	n=len(cpData)
	dists=np.zeros(n)
	for i in range(0,n):
		dists[i]=np.abs(cpData[i]['z']-zPos)
	iout=dists.argmin()
	return cpData[iout]

def getNearestCPData(cpData,zPos):
	n=len(cpData)
	dists=np.zeros(n)
	for i in range(0,n):
		dists[i]=np.abs(cpData[i]['z']-zPos)
	iout=dists.argmin()
	return cpData[iout]

