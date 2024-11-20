import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

testlib = ctypes.CDLL('./testlib.so')
testlib.myprint()
#testlib.getUVW.restype = ctypes.c_void_p
#testlib.getUVW.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double), ctypes.c_double ]

res=np.zeros(3)
pt1=np.array([1,2,3])#np.zeros(3)
pt2=np.array([4,5,6])#np.zeros(3)
pt =np.array([7,8,9])#np.zeros(3).astype(np.float32)

G=ctypes.c_double(10.0)

#testlib.getUVW(np.ascontiguousarray(res, np.float64), np.ascontiguousarray(pt1, np.float64) , np.ascontiguousarray(pt2, np.float64) ,np.ascontiguousarray(pt, np.float64),G)

res=np.array([[0,1,2],[3,4,5],[6,7,8]],dtype=np.float64)


testlib.vortexLine_getUVW.restype = ctypes.c_void_p
testlib.vortexLine_getUVW.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double), ctypes.c_double ]
testlib.vortexRing_getUVW.restype = ctypes.c_void_p
#(double* uvw, double* P1, double* P2 , double* P3 , double* P4 ,double* PT , double G)
testlib.vortexRing_getUVW.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),  ctypes.c_double ]


testlib.testNDArray.restype = ctypes.c_void_p
testlib.testNDArray.argtypes = [ndpointer(dtype=np.float64,ndim=2,shape=res.shape, flags='C_CONTIGUOUS'),ctypes.c_int]

testlib.testNDArray(res,3)

############################# point in polygon algorithm #############################
testlib.vect_pointInPolyXY.restype = ctypes.c_int;
testlib.vect_pointInPolyXY.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ctypes.c_int,ndpointer(ctypes.c_double)]

pts=np.array([ [0,0,0],[1,0,0],[1,1,0],[0,1,0] ],dtype=np.float64)
PT=np.array([0.5,0.5,0.0],dtype=np.float64)

isInPoly=testlib.vect_pointInPolyXY(pts,4,PT)


############################ make CSYS from 2 vectors ############################

testlib.vect_makeCSYS_from2Vectors.restype = ctypes.c_void_p;
testlib.vect_makeCSYS_from2Vectors.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double)]


csys0=np.zeros((3,3),dtype=np.float64)
v1=np.array([1,0,0],dtype=np.float64)
v2=np.array([0,1,0],dtype=np.float64)

testlib.vect_makeCSYS_from2Vectors(csys0,v1,v2)

############################ get transformation matrix #############################

testlib.vect_getTransformMatrix.restype = ctypes.c_void_p;
testlib.vect_getTransformMatrix.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS')]

trMatrix=np.zeros((3,3),dtype=np.float64)
csys2=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]])
testlib.vect_getTransformMatrix(trMatrix,csys2,csys0)

v=np.array((1,0,0))
vTransformed=np.dot(trMatrix,v)


############################ mutliply Two matrices ###################################
#mmath_matrixMultiply(double*res,double*matA,double*matB,int M,int N, int K)
testlib.mmath_matrixMultiply.restype = ctypes.c_void_p;
testlib.mmath_matrixMultiply.argtypes = [ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ctypes.c_int,ctypes.c_int,ctypes.c_int]


#mat1=np.array([[1,2,3],[4,5,6,],[7,8,9]],dtype=np.float64);
#mat2=np.array([1,1,1],dtype=np.float64);
#matmul=np.zeros((3,1),dtype=np.float64)

#testlib.mmath_matrixMultiply(matmul,mat1,np.array([mat2]),mat1.shape[0],mat1.shape[1],1)

############################ transform vector #######################################
#void vect_transformVector(double*res,double*v,double*targetCSYS,double*sourceCSYS){

testlib.vect_transformVector.restype = ctypes.c_void_p;

testlib.vect_transformVector.argtypes = [ndpointer(ctypes.c_double), ndpointer(ctypes.c_double), ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS')]

csys0=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]])
csys2=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]])
trVect=np.array([0,0,0],dtype=np.float64);

v=np.array([1,1,0],dtype=np.float64)

testlib.vect_transformVector(trVect,v,csys2,csys0)

############################ source edge ###############################################
#void sourceEdge_getUVW_local(double *res, double *P1, double* P2, double *PK, double S)
#testlib.sourceEdge_getUVW_local.restype= ctypes.c_void_p;

#testlib.sourceEdge_getUVW_local.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ctypes.c_double]

uvwS=np.zeros(3,dtype=np.float64)

#testlib.sourceEdge_getUVW_local(uvwS,np.array([0,0,0],dtype=np.float64),np.array([10,0,0],dtype=np.float64),np.array([5,5,0],dtype=np.float64),1.)

############################ doubletSourcePanel getUVW ################################################################
#void doubletSourcePanel_getUVW(double*res,double*pts,int n, double*ORIG, double*csysLoc, double*PT,double G,double S)

testlib.sourcePanel_getUVW.restype = ctypes.c_void_p;
testlib.sourcePanel_getUVW.argtypes = [ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ctypes.c_int,
		ndpointer(ctypes.c_double),
		ndpointer(dtype=np.float64,ndim=2, flags='C_CONTIGUOUS'),
		ndpointer(ctypes.c_double),
		ctypes.c_double]
		#ctypes.c_double]
uvw=np.array([10,0,0],dtype=np.float64)
orig=np.array([0.5,0.5,0.],dtype=np.float64)
G=1.0
S=0.0
csys2=np.array([[ 1.,  0.,  0.], [ 0.,  1.,  0.],  [ 0.,  0.,  1.]])
#testlib.sourcePanel_getUVW(uvw,pts,pts.shape[0],orig,csys2,PT,1.0)
