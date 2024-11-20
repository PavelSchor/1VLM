import numpy as np
from numpy import linalg as LA
import math
class B3vect(object):
# pavel.schor@email.cz 2011
#

	def makeVect(self,pt1,pt2):
		if type(pt1) == list or type(pt2) == list:
			pt1=np.array(pt1)
			pt2=np.array(pt2)
		return pt2-pt1

	def vectorMagn(self,v):
		#return np.linalg.norm(v)
	    return math.sqrt(v[0]**2 + v[1]**2 +v[2]** 2)

	def cross(self,left,right):
		return np.array([  	((left[1] * right[2]) - (left[2] * right[1])),
					((left[2] * right[0]) - (left[0] * right[2])),
					((left[0] * right[1]) - (left[1] * right[0]))    ])
	def unitVect(self,v):
		return (v)/np.linalg.norm(v)

	def normalVect(self,u,v):
		return np.cross(u,v)

	def normal1Vect(self,u,v):
		r=np.cross(u,v)
		return r/np.linalg.norm(r)

	def setVectSize(self,v,size):
		vnor=v/self.vectorMagn(v)
		res=vnor*size

		return res

	def getMomentAtPoint(self,p1,p2,F):
		v=self.makeVect(p1,p2)
		M=np.array([0.,0.,0.])
		F1=F
		M[0]= (( F1[1]+F[1])/2.0*v[2]) - ((F1[2]+F[2])/2.0*v[1])
		M[1]=-(( F1[0]+F[0])/2.0*v[2]) + ((F1[2]+F[2])/2.0*v[0])
		M[2]= (( F1[0]+F[0])/2.0*v[1]) - ((F1[1]+F[1])/2.0*v[0])
		return M

	def angle2Vectors(self,v1,v2):
		return (np.dot(v1,v2)/(np.linalg.norm(v1) *np.linalg.norm(v2)))

	def proj(self,u,v):
		return u*np.dot(v,u)/np.dot(u,u)

	def grammSmithOrtho(V):
		U=np.copy(V)
		for i in range(1,V.shape[0]):
			for j in range(0,i):
				U[i] -= self.proj(U[j],V[i])
			U[i] /= np.sqrt( (U[i]**2).sum(axis=0))
		return U

	def makeCSYS2Vectors(self,u,v):
		r=np.zeros((3,3))
		r[0]=u/np.linalg.norm(u)
		r[1]=u/np.linalg.norm(v)
		r[2]=np.cross(r[0],r[1])
		return r
