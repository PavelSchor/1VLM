import numpy as np
import matplotlib.path as mplPath
from utils.loftedSurface import *


#aRot=np.deg2rad(20)
#reverse=False
#pts=wingGeometry.sections[0].sectionA
#plt.plot(pts[:,0],pts[:,2])
#axis=np.array((0,1,0))
#xRudder=0.7
#orig=np.array((xRudder,0,0))
#pts=self.ips[i]
#orig=self.origsGeometry[i]
#axis=self.axis
#aRot=self.deflection
#xRudder=self.teGeometry[i][0]
#reverse=False
def getRudderPtsMask(points,orig,axis,aRot,xRudder,reverse=False):
	pts=np.copy(points)
	p=Plane3D()
	#p.n=np.array((0,1,0))

	p.setFromThreepoints(pts[0], pts[len(pts)/4], pts[len(pts)/2] )
	p.o=np.array([0,0,0])

	pRot=pts[pts[:,0]>xRudder]
	pFix=pts[pts[:,0]<=xRudder]
	rudderPt=pts[:,0]>xRudder
	pt0=p.projectPTS2D(pts)
	bbPath = mplPath.Path(pt0)



	activePts=np.ones(len(pts),dtype=np.bool)
	goesInside=np.zeros(len(pts),dtype=np.bool)

	if not reverse:
		if aRot >=0.:
			for i in range(len(pts)-1,-1,-1):
				if rudderPt[i] == True:
					goesInside[i]=True
				else:
					break
		else:
			for i in range(0,len(pts)):
				if rudderPt[i] == True:
					goesInside[i]=True
				else:
					break

	if reverse:
		if aRot >=0.:
			for i in range(0,len(pts)):
				if rudderPt[i] == True:
					goesInside[i]=True
				else:
					break
		else:
			for i in range(len(pts)-1,-1,-1):
				if rudderPt[i] == True:
					goesInside[i]=True
				else:
					break


	pts[pts[:,0]>xRudder]=bc.rotAx(pRot,axis,orig,aRot)


	pt=p.projectPTS2D(pts)



	ptPlane=p.projectPTS2D(pts)
	ptFixed=p.projectPTS2D(pFix)




	for i in range(0,len(pts)):
		if rudderPt[i]:
			if bbPath.contains_point(ptPlane[i]) and goesInside[i]:
				activePts[i]=False
	return activePts

def getGeometryFromPoints(inFiles,nIterps,mirrorCoord):
	
	geometry=MultiSectionSurface()
	for i in range(0,len(inFiles)-1):
		s0=inFiles[i]

		if type(mirrorCoord) == int:
			s0[:,mirrorCoord]*=-1.

		if type(mirrorCoord) == list:
			for m in mirrorCoord:
				s0[:,m]*=-1.


		s1=inFiles[i+1]


		if type(mirrorCoord) == int:
			s1[:,mirrorCoord]*=-1.

		if type(mirrorCoord) == list:
			for m in mirrorCoord:
				s1[:,m]*=-1.


		geometry.addSection(s0,s1)
		geometry.sections[len(geometry.sections)-1].setNumberOfIterpolants(nIterps[i])
	return geometry

class Rudder(object):
	def __init__(self):
		self.deflection=0.0
		self.axis=np.array((0,1,0))
	
	def setOrigs(self,pts,nIterps,mirrorCoord=None):
		self.origsGeometry= getGeometryFromPoints(pts,nIterps,mirrorCoord).getRawInterpolants()
		
	def setTe(self,pts,nIterps,mirrorCoord=None):
		self.teGeometry= getGeometryFromPoints(pts,nIterps,mirrorCoord).getRawInterpolants()
		
	def setActivePtsMasks(self,reverse=False):
		self.ips=self.geometry.getRawInterpolants()
		self.activePts=np.ones(self.ips[0].shape[0],dtype=np.bool)
		for i in range(0,len(self.ips)):
			self.activePts*=getRudderPtsMask(self.ips[i],self.origsGeometry[i],self.axis,self.deflection,self.teGeometry[i][0],reverse=False)
	
	def deflect(self,res,reverse=False):
		if self.deflection !=0.0:
			self.setActivePtsMasks(reverse)
			
			for i in range(0,len(res)):
					res[i][res[i][:,0]>self.teGeometry[i][0]]=bc.rotAx(res[i][res[i][:,0]>self.teGeometry[i][0]],self.axis,self.origsGeometry[i] ,self.deflection)
		else:
			self.ips=self.geometry.getRawInterpolants()
			self.activePts=np.ones(self.ips[0].shape[0],dtype=np.bool)
