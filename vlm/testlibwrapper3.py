import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from plotPHI import *

testlib = ctypes.CDLL('./testlib.so')

#void vect_transformVector(double*res,double*v,double*targetCSYS,double*sourceCSYS)
testlib.vect_transformVector.restype = ctypes.c_void_p;
testlib.vect_transformVector.argtypes = [
		ndpointer(ctypes.c_double),
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS')]

#void vect_makeCSYS_from2Vectors12(double*res,double*u,double*v){
testlib.vect_makeCSYS_from2Vectors12.restype = ctypes.c_void_p;
testlib.vect_makeCSYS_from2Vectors12.argtypes = [
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ndpointer(ctypes.c_double)]


#void doubletSourcePanel_getUVW(double*res,double*pts,int n, double*ORIG, double*csysLoc, double*PT,double G,double S)
testlib.doubletSourcePanel_getPHI.restype = ctypes.c_double;
testlib.doubletSourcePanel_getPHI.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),	ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_double,
		ctypes.c_double]

#void doubletSourcePanel_getUVW(double*res,double*pts,int n, double*ORIG, double*csysLoc, double*PT,double G,double S)
testlib.doubletSourcePanel_getUVW.restype = ctypes.c_void_p;
testlib.doubletSourcePanel_getUVW.argtypes = [ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_double,
		ctypes.c_double]

#double VXJE_source_edge_influence( double *node_a, double *node_b,double*x)
testlib.VXJE_source_edge_influence.restype = ctypes.c_double;
testlib.VXJE_source_edge_influence.argtypes = [
		ndpointer(ctypes.c_double),
		ndpointer(ctypes.c_double),
		ndpointer(ctypes.c_double)]


testlib.vect_makeGramSmithCS.restype = ctypes.c_void_p;
testlib.vect_makeGramSmithCS.argtypes = 	[ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
						ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS')]

uvw=np.array([0,0,0],dtype=np.float64)
orig=np.array([0.5,0.5,0.],dtype=np.float64)
G=1.0
S=0.0
csys=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]])
pts=np.array([ [0,0,0],[1,0,0],[1,1,0],[0,1,0] ],dtype=np.float64)

PT=np.array([0.5,0.5,1.0000])
#orig=np.array([  9.00000000e-01,  -1.09350000e-02,  -4.16666667e+03])


phip= testlib.doubletSourcePanel_getPHI(pts,pts.shape[0],orig,csys,PT.astype(np.float64),1.0,0.0)

print "phi= ", phip

testlib.doubletSourcePanel_getUVW(uvw,pts,pts.shape[0],orig,csys,PT.astype(np.float64),1.0,0.0)
print uvw;
##0.000000 -0.050240 0.050240 
##1.666667 1.666667 1.666667 
##0.000000 -0.000020 -0.000020 
#pts=np.array([ [-0.050240,1.666667,-0.000020],[0.050240,1.666667,-0.000020],[1,1.666667,-0.000020],[0,1.666667,-0.000020] ],dtype=np.float64)

#PT=np.array([0.,1.6666667,0.000])

#plotZeroPSI_plane( testlib.doubletSourcePanel_getPHI, pts, orig, csys ,G=1.0,S=0.0,plane='xz',c=0.0)
#plotVel_plane( testlib.doubletSourcePanel_getUVW, pts, orig, csys ,G=3.0,S=0.0,plane='xz',c=0.0)



csys1=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]],dtype=np.float64)
csys2=np.array([[ 1.,  0.,  0.], [ 0.,  -1.,  0.],  [ 0.,  0.,  -1.]],dtype=np.float64)
vr=np.array([0,0,0],dtype=np.float64)
v=np.array([0.70710678118654757,0.70710678118654757,0.],dtype=np.float64)

testlib.vect_transformVector(vr,v,csys2,csys1)



csys3=np.zeros((3,3),dtype=np.float64)
vx=np.array([1,0,0],dtype=np.float64)
vy=np.array([0,1,0],dtype=np.float64)

testlib.vect_makeCSYS_from2Vectors12(csys3,vx,vy)


nodeA=np.array([ -0.113114, -0.500000, 0.000000 ])
nodeB=np.array([ -0.113114,  0.500000 ,0.000000 ])
X=np.array([ -0.263592, 0.000000, 0.068058 ])
infl=testlib.VXJE_source_edge_influence(nodeA,nodeB,X)

c=np.cos(np.pi/4.);s=np.sin(np.pi/4.);v=np.array([ [c,s,0],[0,1,0],[0,0,1] ]);

csys4=np.zeros((3,3),dtype=np.float64);
testlib.vect_makeGramSmithCS(csys4,v);


