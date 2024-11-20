import numpy as np
#from scipy.optimize import curve_fit
import time
import pylab as plt
from utils.bCad import *
from scipy.spatial import ConvexHull
#from panelMethod.panelMethod1 import *
#from utils.beamCalc import BeamCrossSection
#from fem.fem4 import FEM1
#from fem.beam3d_nl import Beam3D


def pitchingMomentCoeffFunction(fargs):
	PG=fargs['PG']
	alpha=fargs['alpha']
	q=fargs['q']
	S=fargs['S']
	c=fargs['cmac']
	f=PG.getForceAtPT(fargs['refPT'])[4]
	return f/q/S/c


def liftFunction(f,alpha):
	return -f[0]*np.sin(alpha)+f[2]*np.cos(alpha)
def dragFunction(f,alpha):
	return f[0]*np.cos(alpha)+f[2]*np.sin(alpha)

def liftCoeffFunction(fargs):
	PG=fargs['PG']
	alpha=fargs['alpha']
	q=fargs['q']
	S=fargs['S']
	ix=fargs['ix']
	iz=fargs['iz']
	f=PG.getForceAtRefPT()
	return (-f[0]*np.sin(alpha)+f[2]*np.cos(alpha))/q/S

def dragCoeffFunction(fargs):
	PG=fargs['PG']
	alpha=fargs['alpha']
	q=fargs['q']
	S=fargs['S']
	ix=fargs['ix']
	iz=fargs['iz']
	f=PG.getForceAtRefPT()
	return (f[0]*np.cos(alpha)+f[2]*np.sin(alpha))/q/S

def pitchingMomentCoeffFunction(fargs):
	PG=fargs['PG']
	alpha=fargs['alpha']
	q=fargs['q']
	S=fargs['S']
	c=fargs['cmac']
	f=PG.getForceAtPT(fargs['refPT'])[4]
	return f/q/S/c


#n=section.getNumberOfLiftPanels()
#npan=int(n)
##area=np.zeros(nsec)
#pts=np.zeros((npan*4,3))
#for i in range(0,npan):
	#pan=section.getLiftPanel(i)
	#pts[i*4:i*4+4 , :]=np.copy(pan.pts)
	#if i==1:
		#break
###################################################

def examineSubsectionPointsPosition(section,nsec):
	#n=section.getNumberOfLiftPanels()/float(nsec)
	#if abs(float(n)-float(int(n)))>1e-6:
		#print('Warning number of points does not match with number sections!')
		#return 0
	#npan=int(n)
	nsec=section.nSpanwise
	npan=section.nLPanChordwise

	pts=np.zeros((nsec,3))
	f=np.zeros((npan,3))
	for s in range(0,nsec):
		f.fill(0.0)
		for i in range(0,npan):
			f[i]=section.getLiftPanel(s*npan + i).getCpoint()
		pts[s][0]=np.median(f[:,0]);
		pts[s][1]=np.median(f[:,1]);
		pts[s][2]=np.median(f[:,2]);
	return pts

def findNSEC1(section,maxSec=1000):
	ntry=0
	i=3
	while i<maxSec and ntry==0:
		i+=1
		n=section.getNumberOfLiftPanels()/float(i)
		if abs(float(n)-float(int(n)))<1e-11:
			ntry=i
	return ntry

def findNSEC(section,maxSec=1000):
	ntry=[]
	i=3
	for nsec in range(3,maxSec):
		n=section.getNumberOfLiftPanels()/float(nsec)
		if abs(float(n)-float(int(n)))<1e-11:
			ntry.append(nsec)
	return ntry


def getConvexHullSections(section,pl):
	#nsec=section.nSpanwise
	#n=section.getNumberOfLiftPanels()/float(nsec)
	#if abs(float(n)-float(int(n)))>1e-6:
		#print('Warning number of points does not match with number sections!')
		#return 0
	#npan=int(n)
	
	nsec=section.nSpanwise
	npan=section.nLPanChordwise

	area=np.zeros(nsec)
	ch={}
	for s in range(0,nsec):
		pts=np.zeros((npan*4,3))
		for i in range(0,npan):
			pan=section.getLiftPanel(s*npan + i)
			pts[i*4:i*4+4]=np.copy(pan.pts)
		ch[s]=convexHull2DPoints(pts,pl)
	return ch

def examineSubsectionAreas(section,pl):
	nsec=section.nSpanwise
	ch=getConvexHullSections(section,pl)
	area=np.zeros(nsec)
	for i in ch:
		area[i]=polyArea(ch[i])
	return area


def examineSubsectionForce(section,pts):
	#nsec=len(pts)
	#n=section.getNumberOfLiftPanels()/float(nsec)
	#if abs(float(n)-float(int(n)))>1e-6:
		#print('Warning number of points does not match with number sections!')
		#return 0
	#npan=int(n)
	
	nsec=section.nSpanwise
	npan=section.nLPanChordwise

	cf=np.zeros((nsec,3))
	f=np.zeros(3)
	for s in range(0,nsec):
		f.fill(0.0)
		for i in range(0,npan):
			f+=section.getLiftPanel(s*npan + i).getForce()
		cf[s]=f
	return cf

def examineSubsectionMoment(section,pts):
	#nsec=len(pts)
	#n=section.getNumberOfLiftPanels()/float(nsec)
	#if abs(float(n)-float(int(n)))>1e-6:
		#print('Warning number of points does not match with number sections!')
		#return 0
	#npan=int(n)
	
	nsec=section.nSpanwise
	npan=section.nLPanChordwise

	cm=np.zeros((nsec,3))
	m=np.zeros(3)
	for s in range(0,nsec):
		m.fill(0.0)
		for i in range(0,npan):
			p=section.getLiftPanel(s*npan + i)
			r=p.getCpoint()-pts[s]
			m+=np.cross(r,p.getForce())
		cm[s]=m
	return cm

def examineSubsectionLift(section,pts,liftFunction,alpha,q,areas):
	nsec=section.nSpanwise
	npan=section.nLPanChordwise
	cl=np.zeros(nsec)
	f=np.zeros(3)
	for s in range(0,nsec):
		f.fill(0.0)
		for i in range(0,npan):
			f+=section.getLiftPanel(s*npan + i).getForce()
		cl[s]=liftFunction(f,alpha)#/q/areas[s]
	return cl

def findSubsectionSparPoints(section,spar,j=1,oi=10,i1=50,i2=80):
	nsec=section.nSpanwise
	npan=section.nLPanChordwise
	pts=np.zeros((nsec,3))
	for s in range(0,nsec):
		o=section.getLiftPanel(s*npan + oi).pts[j]
		p1=section.getLiftPanel(s*npan + i1).pts[j]
		p2=section.getLiftPanel(s*npan + i2).pts[j]
		pll=Plane3D()
		pll.setFromThreepoints(o,p1,p2)
		pts[s]=pll.getPolyLineIntersect(spar)
	return pts[np.all(np.isfinite(pts),axis=1)]

def extractSpanwisePanelGroupsFromSection(section):
	nsec=section.nSpanwise
	npan=section.nPanChordwise
	pg={}
	for s in range(0,nsec):
		pg[s]=PanelGroup()
		for i in range(0,npan):
			pg[s].addPanel(section.getPanel(s*npan + i))
	return pg

#def extractCutFromSingleSection(section,j1,j2,includeWake=False):
	

def findSubsectionMassDistributionByAreas(section,pl,mass):
	nsec=section.nSpanwise
	npan=section.nLPanChordwise

	area=np.zeros(nsec)
	cg={}
	pts=np.zeros((npan*4,3))
	for s in range(0,nsec):
		pts.fill(0.0)
		for i in range(0,npan):
			pan=section.getLiftPanel(s*npan + i)
			pts[i*4:i*4+4]=np.copy(pan.pts)
		cg[s]=np.array([pts[:,0].mean(),pts[:,1].mean(),pts[:,2].mean()])
	area=examineSubsectionAreas(section,pl)
	return cg,area/np.sum(area)*mass


def examineSubsectionLDMAtSpar(section,spar,pl,liftFunction,dragFunction,alpha,q,j=1,oi=10,i1=50,i2=80,userSpar=None):
	pts=findSubsectionSparPoints(section,spar,j,oi,i1,i2)
	if userSpar is not None:
		pts=userSpar
	areas=examineSubsectionAreas(section,pl)
	chords=np.zeros(len(pts))
	for i in range(1,len(pts)):
		b=LA.norm(pts[i]-pts[i-1])
		chords[i]=areas[i]/b
	chords[0]=areas[i]/LA.norm(pts[1]-pts[0])
	r=np.zeros((len(pts),10))
	r[:,3]=areas
	r[:,4]=chords	
	r[:,5]=examineSubsectionLift(section,pts,liftFunction,alpha,q,areas)/(q*areas)
	r[:,6]=examineSubsectionLift(section,pts,dragFunction,alpha,q,areas)/(q*areas)
	for i in range(0,3):
		r[:,i]=pts[:,i]
		r[:,i+7]=examineSubsectionMoment(section,pts)[:,i]/(q*areas*chords)

	return r
	
def plotSubsectionLift(section,pts,liftFunction,alpha,q,areas,name,xi):
	cl=examineSubsectionLift(section,pts,liftFunction,alpha,q,areas)
	plt.clf()
	plt.plot(xi,cl)
	plt.grid(True)
	plt.savefig(name)
	plt.close()

def extractSubsections(section,pts):
	#nsec=len(pts)
	#n=section.getNumberOfLiftPanels()/float(nsec)
	#if abs(float(n)-float(int(n)))>1e-6:
		#print('Warning number of points does not match with number sections!')
		#return 0
	#npan=int(n)
	nsec=section.nSpanwise
	npan=section.nLPanChordwise

	sec={}
	for s in range(0,nsec):
		sec[s]=subSection()
		for i in range(0,npan):
			sec[s].addPanel(section.getLiftPanel(s*npan + i))
	return sec

def extractSubsectionCP(s,xc=0):
	x=np.zeros(len(s.getLiftPanelList()))
	y=np.zeros(x.shape)
	for i in s.getLiftPanelList():
		p=s.getLiftPanel(i)
		x[i]=p.getCpoint()[xc]
		y[i]=p.getPressureCoef()
	return x,y
		
def line(x,a,b):    
  return a+b*x

def writelog(pathRes):
  localtime = time.asctime( time.localtime(time.time()) )
  f=open(pathRes+'/log.txt','w')
  f.write('# === Panel method result file ===\n')
  f.write('# \tResults generated: '+localtime+'\n\n')

  f.write('# Reference values\n')
  f.write(' \tS='+str(area)+'\n\tCMAC='+str(chord)+'\n\tXMAC='+str(CG0[0])+'\n\tv='+str(sp)+'\n\trho='+str(rho)+'\n')

  f.write('# Aerodynamic center\n')
  f.write(' \tdCM/dCL='+str(dcmdcl)+'\n\tXAC='+str(CG0[0]-chord*dcmdcl)+'\n')

  f.write('\n# Lift data:\n')
  np.savetxt(f,np.hstack((np.rad2deg(alpha_rigid)[np.newaxis].T, CLU_aircraft[np.newaxis].T)))

  #f.write('\n# Drag data:\n')
  #np.savetxt(f,np.hstack((np.rad2deg(alpha_rigid)[np.newaxis].T, CDU_aircraft[np.newaxis].T)))

  f.write('\n# Pithing moment data at CG0:\n')
  np.savetxt(f,np.hstack((np.rad2deg(alpha_rigid)[np.newaxis].T, CMU0_aircraft[np.newaxis].T)))

  f.write('\n# Pithing moment data at CG:\n')
  np.savetxt(f,np.hstack((np.rad2deg(alpha_rigid)[np.newaxis].T, CMU_aircraft[np.newaxis].T)))

  f.write('\n# CLU CMU data at CG0:\n')
  np.savetxt(f,np.hstack((CLU_aircraft[np.newaxis].T, CMU0_aircraft[np.newaxis].T)))

  f.close()
  
def setBeamElementPropsAlongDir(pr1,propBending,propTorsion,c=1,cyy=1.,cxx=1.,ctt=1.):
	interactors=pr1.elastic.interactors
	elem=pr1.elem
	fem=pr1.fem
	properties=pr1.properties
	materials=pr1.materials
	
	for i in interactors:


		li=abs(interactors[i].orig[c])#;print li;
		A=propBending.getValues(li,'Area',kind='linear')
		Ixx=propBending.getValues(li,'Ixx',kind='linear')*cxx
		Iyy=propBending.getValues(li,'Iyy',kind='linear')*cyy
		J=propTorsion.getValues(li,'J',kind='linear')*ctt
		ofsn=propBending.getValues(li,'offsetSparCG')
		properties[i]={ 'elements':[i],'type':'BEAM3D','property':{'materialID':0, 'A':abs(A), 'Iy':abs(Ixx), 'Iz':abs(Iyy) ,'J':abs(J)}   }

		elem[i]=Beam3D()
		elem[i].propertiesID=i
		#setProps(,E,A,G,J,Iy,Iz)
		elem[i].setProps(E=materials[0]['E'],A=abs(A), G=abs(materials[0]['G']), J=abs(J),Iy=abs(Ixx), Iz=abs(Iyy))
		elem[i].offsetNeutralAxis=ofsn
		ids=interactors[i].ids
		(i1,i2)=(ids[0],ids[1])
		elem[i].setNodesByID(fem.NODE,i1,i2,i3=None,oVector=interactors[i].oVector)
		fem.addElement(elem[i],'wing')
		
	fem.materials=materials
	fem.properties=properties
	
  
  

  
def setBeamElementPropsAlongDir0(pr1,l,H,T,B,c=1):
	interactors=pr1.elastic.interactors
	elem=pr1.elem
	fem=pr1.fem
	properties=pr1.properties
	materials=pr1.materials
	bcs= BeamCrossSection()
	
	for i in interactors:


		li=abs(interactors[i].orig[c])
		Hi=np.interp(li,l,H)
		Bi=np.interp(li,l,B)
		Ti=np.interp(li,l,T)
		JK=interactors[i].getReferenceSection().getJk(3.5e-3)
		
		bcs.setBHbh(Bi,Hi,Bi,Hi-(2.*Ti))
		#	print 'HB :', li, Hi, Bi, hi, bi
		#	print 'PROPS:', B.getArea(), B.getIy(), B.getIz(), B.B, B.H, B.b, B.h
		properties[i]={ 'elements':[i],'type':'BEAM3D','property':{'materialID':0, 'A':bcs.getArea(), 'Iy':bcs.getIy(), 'Iz':bcs.getIz() ,'J':JK}   }

		elem[i]=Beam3D()
		elem[i].propertiesID=i
		#setProps(,E,A,G,J,Iy,Iz)
		elem[i].setProps(E=materials[0]['E'],A=bcs.getArea(), G=materials[0]['G'], J=JK,Iy=bcs.getIy(), Iz=bcs.getIz())
		ids=interactors[i].ids
		(i1,i2)=(ids[0],ids[1])
		elem[i].setNodesByID(fem.NODE,i1,i2,i3=None,oVector=interactors[i].oVector)
		fem.addElement(elem[i])
		
	fem.materials=materials
	fem.properties=properties
	
def extractPanelsToNewGroup(PG,wing,rng1,rng2=None,rngSpanwise=None,increaseCounter=False):
	if rngSpanwise==None:
		rngSpanwise=range(0,wing.nSpanwise)
	
	for i in rngSpanwise:
		nChordwise=0
		ofsRow=i*wing.nLPanChordwise
		for j in rng1:
			PG.addPanel(wing.getLiftPanel(j+ofsRow))
			nChordwise+=1
		if rng2!=None:
			for j in rng2:
				PG.addPanel(wing.getLiftPanel(j+ofsRow))
				nChordwise+=1
	if increaseCounter:			
		PG.nSpanwise+=len(rngSpanwise)
	else:
		PG.nSpanwise=len(rngSpanwise)
	PG.nLPanChordwise=nChordwise
	PG.nPanChordwise=nChordwise

def setPressureCoefPG(PG,cp):
	for i in PG.getLiftPanelList():
		p=PG.getLiftPanel(i)
		p.pressureCoef=cp
