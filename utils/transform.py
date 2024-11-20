import numpy as np

class EulerMatrix(object):
	def __init__(self,a,b,c):
		self.a=a
		self.b=b
		self.c=c
	
	def setA(self,a):
		self.a=a
	def setB(self,b):
		self.b=b
	def setC(self,c):
		self.c=c
	def transform(self,VEC):
		sin=np.sin
		cos=np.cos
		a=self.a
		b=self.b
		c=self.c
		sina=sin(a)
		cosa=cos(a)
		sinb=sin(b)
		cosb=cos(b)
		sinc=sin(c)
		cosc=cos(c)
		RM=np.array(\
		[\
		[cosb*cosa, sinc*sinb*cosa - cosc*sina , cosc*sinb*cosa + sinc*sina],\
		[cosb*sina, sinc*sinb*sina + cosc*cosa, cosc*sinb*sina - sinc*cosa],\
		[- sinb, sinc*cosb, cosc*cosb],\
		])
		return np.dot(VEC,RM)


def rotateVectorByAngles(vec,angles):
	sin=np.sin
	cos=np.cos
	a=angles[0]
	b=angles[1]
	c=angles[2]
	sina=sin(a)
	cosa=cos(a)
	sinb=sin(b)
	cosb=cos(b)
	sinc=sin(c)
	cosc=cos(c)
	RM=np.array(\
	[\
	[cosb*cosa, sinc*sinb*cosa - cosc*sina , cosc*sinb*cosa + sinc*sina],\
	[cosb*sina, sinc*sinb*sina + cosc*cosa, cosc*sinb*sina - sinc*cosa],\
	[- sinb, sinc*cosb, cosc*cosb],\
	])
	return np.dot(vec,RM)