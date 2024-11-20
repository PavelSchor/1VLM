import numpy as np
from utils.b3Vect import B3vect
vect=B3vect()

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

testlib = ctypes.CDLL('panelMethod/testlib.so')
testlib.vortexLine_getUVW.restype = ctypes.c_void_p
testlib.vortexLine_getUVW.argtypes = [ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double),ndpointer(ctypes.c_double), ctypes.c_double ]

class vortexPolyLine(dict):

	def __init__(self):
		self.__ns=-1

	def append(self,vl):
		self.__ns+=1
		self[self.__ns]=vl

	def addSegment(self,P1,P2):
		self.append(vortexLine(P1,P2))
	def setG(self,G):
		for i in range(0,len(self)):
			self[i].setG(G)
	def getUVW(self,PT):
		uvw=np.zeros(3)
		for i in range(0,len(self)):
			uvw+=self[i].getUVW(PT)
		return uvw


class vortexLine(object):
	def __init__(self,P1,P2):
		self.G=1.0
		self.P1=P1
		self.P2=P2
		self.PM=self.P1+0.5*(self.P2-self.P1)#vect.makeVect(P1,P2)

	def setG(self,G):
		self.G=float(G)

	def getP1(self):
		return self.P1

	def getP2(self):
		return self.P2

	def getPM(self):
		return self.PM


	def getG(self):
		return self.G


	def getUVW(self,PT):
		err=0.000001
		P1=self.getP1()
		P2=self.getP2()
		PM=self.getPM()
		G=self.getG()
		r0=self.P2-self.P1#vect.makeVect(P1,P2)
		r1=self.P1-PT#vect.makeVect(PT,P1)
		r2=self.P2-PT#vect.makeVect(PT,P2)
		r12=vect.cross(r1,r2)
	
		r1m=vect.vectorMagn(r1)
		r2m=vect.vectorMagn(r2)
		r12m=vect.vectorMagn(r12)
		
		
		uvw=np.zeros(3)
		if r1m> err and r2m > err or r12m > err :
			K=(np.dot(r0,r1)/r1m - np.dot(r0,r2)/r2m)/4.0/np.pi/(r12m*r12m)
			uvw=r12*K
	
			#uvw=r12*(G/4./np.pi/r12m**2*(np.dot(r0,r1)/r1m - np.dot(r0,r2)/r2m) )
		return uvw
#		res=np.zeros(3)

#		testlib.vortexLine_getUVW(np.ascontiguousarray(res, np.float64), np.ascontiguousarray(self.P1, np.float64) , np.ascontiguousarray(self.P2, np.float64) ,np.ascontiguousarray(PT, np.float64),ctypes.c_double(self.G))

#		testlib.vortexLine_getUVW(res, self.P1 , self.P2 ,PT, ctypes.c_double(self.G))

#		return res		
		

