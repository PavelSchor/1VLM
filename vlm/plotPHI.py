import numpy as np
import matplotlib.pyplot as plt
from numpy import arange
from numpy import meshgrid

def plotZeroPSI_plane(fPHI,pts,orig,csys,G=0.0,S=0.0,plane='xz',c=0.0):
	delta = 0.015
	xrg = arange(-2.0, 2.0, delta)
	yrg = arange(-2.0, 2.0, delta)
	Y, X = meshgrid(xrg,yrg)
	F=np.zeros(X.shape)
	H=np.zeros(X.shape)
	U=np.zeros(X.shape)
	V=np.zeros(X.shape)
	for i in range(0,xrg.shape[0]):
		x=xrg[i]
		for j in range(0,yrg.shape[0]):
			y=yrg[j]
			if plane=='xy':
				PT=np.array([x,y,c])
			if plane=='xz':
				PT=np.array([x,c,y])
			if plane=='yz':
				PT=np.array([c,x,y])
			F[i][j] = fPHI(pts,pts.shape[0],orig,csys,PT.astype(np.float64),G,S)
#			U[i][j]=vel[0]
#			V[i][j]=vel[1]

#	U[np.abs(U) > 50] = 0.0
#	V[np.abs(V) > 50] = 0.0
	
	
	#dom.getPHIAtPoint(np.array([x,y,z]))
	#p1=plt.contour(Y, X, (F), np.linspace(0,5,100))
	#p2=plt.contour(Y, X, (F), np.linspace(0,-10,100))
	#p3=plt.contour(Y, X, (F), [0],linewidths=2)
	
	p3=plt.contour(Y, X, (F),  np.linspace(-5,5,200))

#	p3.collections[0].set_label(r'$\psi = 0$')
	
#	plt.quiver(X,Y,U,V,units='xy')
	#plt.quiver(U,V)


	plt.legend(loc='upper left')
	plt.grid(True)
	plt.xlabel(r'$ x$')
	plt.ylabel(r'$ y$')
	plt.show()
	return X,Y,F#U,V

def plotVel_plane(fPHI,pts,orig,csys,G=0.0,S=0.0,plane='xz',c=0.0):
	delta = 0.01
	xrg = arange(-2.0, 2.0, delta)
	yrg = arange(-2.0, 2.0, delta)
	Y, X = meshgrid(xrg,yrg)
	F=np.zeros(X.shape)
	H=np.zeros(X.shape)
	U=np.zeros(X.shape)
	V=np.zeros(X.shape)
	for i in range(0,xrg.shape[0]):
		x=xrg[i]
		for j in range(0,yrg.shape[0]):
			y=yrg[j]
			uvw=np.zeros(3,dtype=np.float64)
			if plane=='xy':
				PT=np.array([x,y,c])
				fPHI(uvw,pts,pts.shape[0],orig,csys,PT.astype(np.float64),G,S)
				U[i][j]=uvw[0]
				V[i][j]=uvw[1]
			if plane=='xz':
				PT=np.array([x,c,y])
				fPHI(uvw,pts,pts.shape[0],orig,csys,PT.astype(np.float64),G,S)
				U[i][j]=uvw[0]
				V[i][j]=uvw[2]
			if plane=='yz':
				PT=np.array([c,x,y])
				fPHI(uvw,pts,pts.shape[0],orig,csys,PT.astype(np.float64),G,S)
				U[i][j]=uvw[1]
				V[i][j]=uvw[2]

#			F[i][j] = fPHI(pts,pts.shape[0],orig,csys,PT.astype(np.float64),G,S)
#			U[i][j]=vel[0]
#			V[i][j]=vel[1]

#	U[np.abs(U) > 50] = 0.0
#	V[np.abs(V) > 50] = 0.0
	
	
	#dom.getPHIAtPoint(np.array([x,y,z]))
	#p1=plt.contour(Y, X, (F), np.linspace(0,5,100))
	#p2=plt.contour(Y, X, (F), np.linspace(0,-10,100))
	#p3=plt.contour(Y, X, (F), [0],linewidths=2)
	
#	p3=plt.contour(Y, X, (F),  np.linspace(-5,5,100))

#	p3.collections[0].set_label(r'$\psi = 0$')
	
	plt.quiver(X,Y,U,V)#,units='dots')
	#plt.quiver(U,V)


	plt.legend(loc='upper left')
	plt.grid(True)
	plt.xlabel(r'$ x$')
	plt.ylabel(r'$ y$')
	plt.show()
	return X,Y,U,V
