from scipy.spatial import ConvexHull
import numpy as np
import math

from utils.b3Vect import *

LA=np.linalg

def commonEdgeTwoPanels(A,B):
	nrows, ncols = A.pts.shape
	dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [A.pts.dtype]}
	C = np.intersect1d(A.pts.view(dtype), B.pts.view(dtype))
	if len(C) < 2:
		return False
	else:
		return True

def findLeftRightCloseToZero(a):
	left=a.min()
	right=a.max()
	li=a.argmin()
	ri=a.argmax()
	for i in range(0,len(a)):
		if a[i]>0:
			if a[i] < right:
				right =a[i]
				ri=i
		if a[i]<0:
			if a[i] > left :
				left =a[i]
				li=i
	return li,ri

class Plane3D(object):
	def __init__(self):
		self.o=np.array((0.,0.,0.))
		self.n=np.array((0.,0.,1.))
		self.x=np.array((1.,0.,0.))
		self.y=np.array((0.,1.,0.))
	
	def setN(self,n):
		self.n=n/LA.norm(n)
	
	def setO(self,o):
		self.o=o
	
	def setFromThreepoints(self,o,p1,p2):
		self.setO(o)
		self.x=p1-o
		self.x/=LA.norm(self.x)
		self.y=p2-o
		self.y/=LA.norm(self.y)
		self.n=np.cross(self.x,self.y)
		self.n/=LA.norm(self.n)
		
	def projectPT(self,pt):
		
		#1) Make a vector from your orig point to the point of interest:
		v = pt-self.o
		#2) Take the dot product of that vector with the unit normal vector n:

		dist=np.dot(v,self.n)
		#3) Multiply the unit normal vector by the distance, and subtract that vector from your point.
		return pt - dist*self.n
	
	def projectPTS2D(self,pts):
		r=np.zeros((len(pts),2))
		for i in range(0,len(pts)):
			 r[i]=self.projectPT2D(pts[i])
		return r
	
	def projectPT2D(self,pt):
		p=self.projectPT(pt)
		x = np.dot((p - self.o), self.x)
		y=np.dot((p - self.o), np.cross(self.n,self.x))
		return np.array((x,y))
	
	def vector2x(self,v):
		x=self.projectPT(v)
		self.x=x/LA.norm(x)

	def vector2y(self,v):
		y=self.projectPT(v)
		self.y=y/LA.norm(y)
		
	def getPointDistance(self,pt):
		#http://paulbourke.net/geometry/pointlineplane/
		#d=(A xa + B ya + C za + D) / sqrt(A2 + B2 + C2)
		#distance = (Pa - Pb) dot n / ||n||
		return np.dot(pt-self.o, self.n)/LA.norm(self.n)

	def getLineIntersect(self,l):
		#http://paulbourke.net/geometry/pointlineplane/
		#P = P1 + u (P2 - P1)
		#u = (n dot(o - p1) / (n dot(p2-p1))
		p1=l[0]
		p2=l[1]
		u=np.dot(self.n,self.o-p1)/np.dot(self.n,p2-p1)
		return p1+u*(p2-p1)

	def getPolyLineIntersect(self,p):
		d=np.zeros(p.shape[0])
		for i in range(0,p.shape[0]):
			d[i]=self.getPointDistance(p[i])
		l,r=findLeftRightCloseToZero(d)
		line=np.array([ p[l],p[r]])
		return self.getLineIntersect(line)
	
	
# Graham Scan - Tom Switzer <thomas.switzer@gmail.com>
TURN_LEFT, TURN_RIGHT, TURN_NONE = (1, -1, 0)

def turn(p, q, r):
    return cmp((q[0] - p[0])*(r[1] - p[1]) - (r[0] - p[0])*(q[1] - p[1]), 0)

def _keep_left(hull, r):
    while len(hull) > 1 and turn(hull[-2], hull[-1], r) != TURN_LEFT:
            hull.pop()
    if not len(hull) or hull[-1] != r:
        hull.append(r)
    return hull

def convex_hull(pts):
    """Returns points on convex hull of an array of points in CCW order."""
    points=pts.tolist()
    points = sorted(points)
    l = reduce(_keep_left, points, [])
    u = reduce(_keep_left, reversed(points), [])
    res= l.extend(u[i] for i in xrange(1, len(u) - 1)) or l
    return np.array(res)

def polyArea1(poly):
	A=0.0
	for i in range(0,poly.shape[0]-1):
		A+= (poly[i,0]*poly[i+1,1] - poly[i+1,0]*poly[i,1])
		return 0.5*A

def polyArea(poly):
	x=poly[:,0]
	y=poly[:,1]
	return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))


def convexHull2DPoints(points,pl):
	pts=np.zeros((points.shape[0],2))
	for i in range(pts.shape[0]):
		pts[i]=pl.projectPT2D(points[i])
	cv= ConvexHull(pts)
	return pts[cv.vertices]

def interpByCoords(i,xmin,xmax,pts,ni=100):
	nr,nc=pts.shape
	if nc == 2:
		inp=np.hstack((pts,np.zeros(nr)))
	else:
		inp=pts
	r=np.zeros((ni,3))
	if abs(xmin) <= abs(xmax):
		x=abs(np.linspace(xmin,xmax,ni))
	else:
		x=abs(np.linspace(xmax,xmin,ni))
	for j in xrange(0,3):
		r[:,j]=np.interp(x,abs(inp[:,i]),inp[:,j])
	return r

class BCad(object):
# 2011 pavel.schor@email.cz
#
	def pointsDist(self,P1,P2):
		return np.sqrt( (P2[0]-P1[0])**2.0 + (P2[1]-P1[1])**2.0 + (P2[2]-P1[2])**2.0 )

	def lineLength(self,l):
		vect=B3vect()
		llen=np.zeros((l.shape[0]-1))
		for i in range(0,len(l)-1 ):
			llen[i]=vect.vectorMagn(vect.makeVect(l[i],l[i+1]) )
		return llen

	def rotX(self,inp,org,fi):
		r=np.zeros((inp.shape[0] ,4))
		org=np.append(org,1)
		for i in range(0, inp.shape[0]):
			pt=inp[i,:]
			pt=np.append(pt,1)
			de=org-pt
			npt=pt-org
			rx=np.array([[1.,0.,0.,0.], [0.,math.cos(fi),math.sin(fi),0.],[0.,-math.sin(fi),math.cos(fi),0],[0,0,0,1]])
			rr=np.dot(npt,rx)+org.transpose() 
			r[i]=rr
		return r[:,0:3]


	def rotY(self,inp,org,fi):
		r=np.zeros((inp.shape[0] ,4))
		org=np.append(org,1)
		for i in range(0, inp.shape[0]):
			pt=inp[i,:]
			pt=np.append(pt,1)
			de=org-pt
			npt=pt-org
			ry=np.array([ [math.cos(fi),0.,-math.sin(fi),0.],[0.,1.,0.,0.],[math.sin(fi),0.,math.cos(fi),0],[0,0,0,1]])
			rr=np.dot(npt,ry)+org.transpose() 
			r[i]=rr
		return r[:,0:3]

	def rotZ(self,inp,org,fi):
		r=np.zeros((inp.shape[0] ,4))
		org=np.append(org,1)
		for i in range(0, inp.shape[0]):
			pt=inp[i,:]
			pt=np.append(pt,1)
			de=org-pt
			npt=pt-org
			rz=np.array([ [math.cos(fi),math.sin(fi),0.,0.],[-math.sin(fi),math.cos(fi),0,0],[0,0,1,0],[0,0,0,1]])
			rr=np.dot(npt,rz)+org.transpose() 
			r[i]=rr
		return r[:,0:3]


	def rotAx(self,inp,ax,org,r):
		#http://answers.google.com/answers/threadview/id/361441.html
		#http://www.cprogramming.com/tutorial/3d/quaternions.html
		u0=ax-org
		u=ax/np.linalg.norm(ax)
		q0 = np.cos(r/2.0); 
		q1 = np.sin(r/2.0)*u[0]
		q2 = np.sin(r/2.0)*u[1]
		q3 = np.sin(r/2.0)*u[2]
		Q=np.array([ 	[q0**2 + q1**2 - q2**2 - q3**2    ,    2.*(q1*q2 - q0*q3)    ,     2.*(q1*q3 + q0*q2) ],
				[2.*(q2*q1 + q0*q3)   ,    (q0**2 - q1**2 + q2**2 - q3**2)   ,     2.*(q2*q3 - q0*q1) ],
				[2.*(q3*q1 - q0*q2)   ,   2.*(q3*q2 + q0*q1)   ,      (q0**2 - q1**2 - q2**2 + q3**2) ] ])
		npt=inp-org
		res=np.zeros((inp.shape[0] ,3))
		for i in range(0, inp.shape[0]):
			pt=npt[i,:]
			res[i]=np.dot(Q,npt[i])+org
		return res
			
	def scl(self,inp,org,s):
		r=np.zeros((inp.shape[0] ,4))
		org=np.append(org,1)
		for i in range(0, inp.shape[0]):
			pt=inp[i,:]
			pt=np.append(pt,1)
			de=org-pt
			npt=pt-org
			sc=np.array([ [s[0],0.,0.,0.],[0.,s[1],0.,0.],[0.,0.,s[2],0.],[0.,0.,0.,1.]])
			rr=np.dot(npt,sc)+org.transpose() 
			r[i]=rr
		return r[:,0:3]

	def transl(self,inp,org,t):
		r=np.zeros((inp.shape[0] ,4))
		org=np.append(org,1)
		for i in range(0, inp.shape[0]):
			pt=inp[i,:]
			pt=np.append(pt,1)
			de=org-pt
			npt=pt-org
			tr=np.array([ [1.,0.,0.,0.],[0.,1.,0.,0.],[0.,0.,1.,0.],[t[0],t[1],t[2],1]])
			rr=np.dot(pt,tr)+org.transpose() 
			r[i]=rr
		return r[:,0:3]

	def intersectByPlane(self,in1,in2,pl):
		r=np.zeros((in1.shape[0] ,3))
		for i in range(0, in1.shape[0]):
			pt1=in1[i,:]
			pt2=in2[i,:]
			A=np.array( [[pt1[0] - pl[0,0]], [ pt1[1] -pl[0,1]], [pt1[2]-pl[0,2]] ])
			B=np.array([ [pt1[0] - pt2[0], pl[1,0]-pl[0,0], pl[2,1]-pl[0,0]],   [pt1[1] - pt2[1], pl[1,1]-pl[0,1], pl[2,1]-pl[0,1]],  [pt1[2] - pt2[2], pl[1,2]-pl[0,2], pl[2,2]-pl[0,2]]   ])
			x=LA.solve(B,A)

			r[i]=pt1+(pt2-pt1)*x[0,0]
		return r


	def lineMidPt(self,l):
		m=l.shape[0]-1
		n=l.shape[1]
		r=np.zeros((m ,n))
		for i in range(0, m):
			for j in range(0,n):
				r[i,j]=(l[i,j] +l[i+1,j])*0.5

		return r


	def ofsXZplane(self,ofs):

		pl=np.array([ [1.,ofs,0.],[2.,ofs,0.],[2.,ofs,1.] ])

		return pl

	def ofsXYplane(self,ofs):

		pl=np.array([ [0.,0.,ofs],[1.,0.,ofs],[1.,1.,ofs] ])

		return pl

	def ofsYZplane(self,ofs):

		pl=np.array([ [ofs,0.,0.],[ofs,1.,0.],[ofs,1.,1.] ])

		return pl


	def polyArea(self,poly,ni=None):
		if ni==None:
			ni=[0,1,2]
		v1=poly[ni[1]]-poly[ni[0]];v2=poly[ni[2]]-poly[ni[0]];
		unitNormal=np.cross(v1,v2)
		unitNormal/=np.linalg.norm(unitNormal)
		if len(poly) < 3:
			return 0
		total = np.zeros(3)
		N = len(poly)
		for i in range(N):
			vi1 = poly[i]
			vi2 = poly[(i+1) % N]
			prod = np.cross(vi1, vi2)
			total += prod
		result = np.dot(total, unitNormal)
		return abs(result/2)

	def polyCircumference(self,poly):
		ln=0.0
		N = len(poly)
		for i in range(N):
		 ln+=np.linalg.norm(poly[(i+1) % N ] - poly[i])
		return ln
	
	def rotatePanels(self,PG,panelList,ax,org,r,mask=np.array([True,True,True,True])):
		for i in panelList:
			p=PG.getPanel(i)
			pts=self.rotAx(p.pts,ax,org,r)
			#p.P1=pts[0]
			#p.P2=pts[1]
			#p.P3=pts[2]
			#p.P4=pts[3]
			for j in range(0,4):
				if mask[j]:
					p.pts[j]=pts[j]
			p.initGeometry()

	def deflectRudder(self,pr1,PG,ax,org,r,pgids,stripMask,spanMask):
		strps=PG.getSpanwiseStrips(True)
		for k in strps:
			s=strps[k]
			for i in pgids:
				p=s.getPanel(i)
				if p.isLast:
					pwo=pr1.dom1.getPanelByGid(p.wakes[0])
					if not hasattr(pwo,'ite1'):
						ptsC=self.commonEdgePointsTwoPanels(p,pwo)
						p1=ptsC[0]
						p2=ptsC[1]
						norm=np.linalg.norm(p.pts-p1,axis=1)
						p.ite1=norm.argsort()[0]
						norm=np.linalg.norm(p.pts-p2,axis=1)
						p.ite2=norm.argsort()[0]


						norm=np.linalg.norm(pwo.pts-p1,axis=1)
						pwo.ite1=norm.argsort()[0]
						pwo.ite3=norm.argsort()[1]

						norm=np.linalg.norm(pwo.pts-p2,axis=1)
						pwo.ite2=norm.argsort()[0]
						pwo.ite4=norm.argsort()[1]


			
			for i in pgids:
				mask=np.array([True,True,True,True])
				if spanMask is not None:
					try:
						mask*=spanMask[k]
					except:
						pass
				try:
					mask*=stripMask[i]
				except:
					pass
				
				p=s.getPanel(i)
				pts=self.rotAx(p.pts,ax,org,r)
				#p.P1=pts[0]
				#p.P2=pts[1]
				#p.P3=pts[2]
				#p.P4=pts[3]
				for j in range(0,4):
					if mask[j]:
						p.pts[j]=pts[j]
				p.initGeometry()
				
			for i in pgids:
				p=s.getPanel(i)
				if p.isLast:
					pwo=pr1.dom1.getPanelByGid(p.wakes[0])
					v1=p.pts[p.ite1]-pwo.pts[pwo.ite1]
					v2=p.pts[p.ite2]-pwo.pts[pwo.ite2]
					for j in p.wakes:
						pw=pr1.dom1.getPanelByGid(j)
						pw.pts[pwo.ite1]+=v1;pw.pts[pwo.ite3]+=v1;
						pw.pts[pwo.ite2]+=v2;pw.pts[pwo.ite4]+=v2;

	def rotatePanelGroup(self,PG,ax,org,r):
		rotatePanels(PG,PG.getPanelList(),ax,org,r)

	def rotatePanelGroupStrip(self,PG,sid,ptl,ax,org,r,liftOnly=True):
		s=PG.getSpanwiseStrips(True)[sid]
		pl=s.getPanelList()
		if not liftOnly:
			try:
				p=PG.getPanel(0)
				pl+=p.wakes
			except:
				pass
		for i in s.getPanelList():
			p=PG.getPanel(i)
			pts=self.rotAx(p.pts,ax,org,r)
			for j in range(0,4):
				if j in ptl:
					p.pts[j]=pts[j]
			p.initGeometry()


	def ptNormal(self,pt,cp1,cp2,cp3,mag):
		vect=B3vect()
		nv={}
		norms={}
		v1=vect.makeVect(cp1,cp3)
		v2=vect.makeVect(cp1,cp2)
		nv=vect.normal1Vect(v1,v2)*mag

		a=np.array([0.,0.,0.])
		b=np.array([0.,0.,0.])

		a[0]=pt[0]
		a[1]=pt[1]
		a[2]=pt[2]

		b[0]=pt[0]+nv[0]
		b[1]=pt[1]+nv[1]
		b[2]=pt[2]+nv[2]

		norms={}
		norms["a"]=a
		norms["b"]=b

		return norms

	def curNormals(self,pt,cp1,cp2,cp3,mag):
		vect=B3vect()
		norms={}
		for i in range(0, cp1.shape[0]):
			cpt1=cp1[i]
			cpt2=cp2[i]
			cpt3=cp3[i]
			ptt=pt[i]
						
			norms[i]= self.ptNormal(ptt,cpt1,cpt2,cpt3,mag)

		return norms


	def transform2Dfoil(self,cpx):
		m=cpx.shape[0]
		n=cpx.shape[1]
		cc=np.zeros( (m,n) )
		for i in range(0, m):

			if cpx[i][0] > cpx[i-1][0]:
				cc[i,:]= cpx[i,:]
				cc[i,0]=-cpx[i][0]
			else:
				cc[i,:]=cpx[i,:]
		cc[0,:]=cpx[0,:]
		
		return cc

	def transform2DfoilBack(self,cpx):
		m=cpx.shape[0]
		n=cpx.shape[1]
		cc=np.zeros( (m,n) )
		for i in range(0, m):

			if cpx[i][0] < 0.:
				cc[i,:]= cpx[i,:]
				cc[i,0]=-cpx[i][0]
			else:
				cc[i,:]=cpx[i,:]
		cc[0,:]=cpx[0,:]
		
		return cc	


	def setFoilCPs(self,pts,cps):
	
		x=pts[0]
		y=pts[1]
		z=pts[2]

		pti=self.lineMidPt(pts)
		ptr=self.transform2Dfoil(pti)

		cptr=self.transform2Dfoil(cps)

		pti=np.sort(pti,axis=0)
		ptr=np.sort(ptr,axis=0)
		cptr=np.sort(cptr,axis=0)
		cpi=np.interp(pti[:,0], cptr[:,0],cptr[:,1])

		
		return cpi


	def getChordLine(self,cut,n):
		
		imin=np.argmin(cut[:,0])
		imax=np.argmax(cut[:,0])
	
		lp=cut[imin,:]
		tp=cut[imax,:]
	
		xo=np.array([lp[0],tp[0]])
		yo=np.array([lp[1],tp[1]])
		zo=np.array([lp[2],tp[2]])
	
		xi=np.linspace(lp[0],tp[0],n)
		
	
		yi=np.interp(xi,xo,yo)

		zi=np.interp(xi,xo,zo)

		
		res=res=np.zeros((n,3))
		
		res[:,0]=xi
		res[:,1]=yi
		res[:,2]=zi
	
		return res



	def findNearestSection(self,direction,PTS,*SECS):
		dmax=np.array(())
		dmin=np.array(())


		for sec in range(0,len(SECS)-1):
			for pt in PTS:
				for spt in sec:
					v=0
					vect=B3vect()
					nv={}
					norms={}
					v1=vect.makeVect(cp1,cp3)
					v2=vect.makeVect(cp1,cp2)
					nv=vect.normal1Vect(v1,v2)*mag

	def commonEdgeTwoPanels(self,A,B):
		nrows, ncols = A.pts.shape
		dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [A.pts.dtype]}
		C = np.intersect1d(A.pts.view(dtype), B.pts.view(dtype))
		if len(C) < 2:
			return False
		else:
			return True

	def commonEdgePointsTwoPanels(self,A,B):
		nrows, ncols = A.pts.shape
		dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [A.pts.dtype]}
		C = np.intersect1d(A.pts.view(dtype), B.pts.view(dtype))
		if self.commonEdgeTwoPanels(A,B):
			r=np.zeros((2,3))
			for i in range(0,3):
				r[0,i]=C[0][i]
				r[1,i]=C[1][i]
			return r
			#return C
		else:
			return None
		
def rotateX(inp,org,fi):
	r=np.zeros((inp.shape[0] ,4))
	org=np.append(org,1)
	for i in range(0, inp.shape[0]):
		pt=inp[i,:]
		pt=np.append(pt,1)
		de=org-pt
		npt=pt-org
		rx=np.array([[1.,0.,0.,0.], [0.,math.cos(fi),math.sin(fi),0.],[0.,-math.sin(fi),math.cos(fi),0],[0,0,0,1]])
		rr=np.dot(npt,rx)+org.transpose() 
		r[i]=rr
	return r[:,0:3]


def rotateY(inp,org,fi):
	r=np.zeros((inp.shape[0] ,4))
	org=np.append(org,1)
	for i in range(0, inp.shape[0]):
		pt=inp[i,:]
		pt=np.append(pt,1)
		de=org-pt
		npt=pt-org
		ry=np.array([ [math.cos(fi),0.,-math.sin(fi),0.],[0.,1.,0.,0.],[math.sin(fi),0.,math.cos(fi),0],[0,0,0,1]])
		rr=np.dot(npt,ry)+org.transpose() 
		r[i]=rr
	return r[:,0:3]

def rotateZ(inp,org,fi):
	r=np.zeros((inp.shape[0] ,4))
	org=np.append(org,1)
	for i in range(0, inp.shape[0]):
		pt=inp[i,:]
		pt=np.append(pt,1)
		de=org-pt
		npt=pt-org
		rz=np.array([ [math.cos(fi),math.sin(fi),0.,0.],[-math.sin(fi),math.cos(fi),0,0],[0,0,1,0],[0,0,0,1]])
		rr=np.dot(npt,rz)+org.transpose() 
		r[i]=rr
	return r[:,0:3]
