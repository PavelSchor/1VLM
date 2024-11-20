import numpy as np
from numpy import linalg as LA
class B3vect(object):
# pavel.schor@email.cz 2011
#

	def makeVect(self,pt1,pt2):
		import numpy as np
		r=np.array([0.,0.,0.])
		r[0]=pt2[0]-pt1[0]
		r[1]=pt2[1]-pt1[1]
		r[2]=pt2[2]-pt1[2]
		
		return r

	def vectorMagn(self,v):
		mag=( v[0]**2. +v[1]**2. + v[2]**2.)**0.5

		return mag

	def unitVect(self,v):
		return (v)/self.vectorMagn(v)

	def normalVect(self,u,v):
		n=np.array([0.,0.,0.])
		n[0] = u[1]*v[2] - u[2]*v[1]
		n[1] = u[2]*v[0] - u[0]*v[2]
		n[2] = u[0]*v[1] - u[1]*v[0]

		return n

	def normal1Vect(self,u,v):
		n=np.array([0.,0.,0.])
		n[0] = u[1]*v[2] - u[2]*v[1]
		n[1] = u[2]*v[0] - u[0]*v[2]
		n[2] = u[0]*v[1] - u[1]*v[0]

		return n/self.vectorMagn(n)


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
