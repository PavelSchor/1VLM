import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

testlib = ctypes.CDLL('./testlib.so')


############################ doubletSourcePanel getUVW ################################################################
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
uvw=np.array([0,0,0],dtype=np.float64)
orig=np.array([0.5,0.5,0.],dtype=np.float64)
G=1.0
S=0.0
csys=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]],dtype=np.float64)
pts=np.array([ [0,0,0],[1,0,0],[1,1,0],[0,1,0] ],dtype=np.float64)
PT=np.array([0.0,0.0,0.0])

testlib.doubletSourcePanel_getUVW(uvw,pts,pts.shape[0],orig,csys,PT.astype(np.float64),1.0,0.0)
print uvw;
