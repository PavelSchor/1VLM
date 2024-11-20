import numpy as np
np.set_printoptions(linewidth=132,precision=9)
import pylab as plt
import copy
import cProfile
from utils.mesh import *
from panelMethod.panelMethod1 import *
from panelMethod.env import GlobalEnvironment
from utils.b3Vect import B3vect 
from utils.bCad import BCad
from resampleFoil import resample

from fem.fem4 import FEM1
from fem.beam3d import Beam3D
from utils.interactor import PanelBeamInteractor
from utils.loftedSurface import *
from utils.elasticWing import ElasticWing1
import utils.csys as cs
from utils.beamCalc import BeamCrossSection as BeamCalc
from fem.utils.writeBDF import BDFWriter
from utils.simpleSTEPReader import extractPointsFromSTEPFile

bdfWriter=BDFWriter()
vect=B3vect()
bc=BCad()


def plotSection(s0,caption=''):
	for i in range(0,len(s0)):
		plt.text(s0[i,1],s0[i,2],str(i))
	plt.plot(s0[:,1],s0[:,2])
	plt.title(caption)
	plt.show()

def makeMultiSectionSurfaceFromArrays(arrays,nIterps,reverse,isReversed,resampling,scale,chord_i=0,mirrorCoord=None,transl=None):
	geometry=MultiSectionSurface()
	for i in range(0,len(arrays)-1):
		s0=arrays[i]

		if type(mirrorCoord) == int:
			s0[:,mirrorCoord]*=-1.

		if type(mirrorCoord) == list:
			for m in mirrorCoord:
				s0[:,m]*=-1.

		if reverse[i]:
			s0=s0[::-1]
		if resampling == 0:
			s0*=scale
		else:
			s0=resample(s0,resampling,chord_i=chord_i,reversed=isReversed[i])*scale

		if transl != None:
			s0=s0+transl

		s1=arrays[i+1]


		if type(mirrorCoord) == int:
			s1[:,mirrorCoord]*=-1.

		if type(mirrorCoord) == list:
			for m in mirrorCoord:
				s1[:,m]*=-1.

		if reverse[i+1]:
			s1=s1[::-1]
		if resampling == 0:
			s1*=scale
		else:
			s1=resample(s1,resampling,chord_i=chord_i,reversed=isReversed[i])*scale

		if transl != None:
			s1+=transl

		geometry.addSection(s0,s1)
		geometry.sections[len(geometry.sections)-1].setNumberOfIterpolants(nIterps[i])
	nSec=len(geometry.getInterpolants())/2
	return geometry,nSec



def makeMultiSectionSurfaceFromSTPFiles(path,inFiles,nIterps,reverse,isReversed,resampling,scale,chord_i=0,mirrorCoord=None,transl=None):
	geometry=MultiSectionSurface()
	for i in range(0,len(inFiles)-1):
		print inFiles[i]
		s0=extractPointsFromSTEPFile(path+inFiles[i])#[1:-1,:]

		if type(mirrorCoord) == int:
			s0[:,mirrorCoord]*=-1.

		if type(mirrorCoord) == list:
			for m in mirrorCoord:
				s0[:,m]*=-1.

		if reverse[i]:
			s0=s0[::-1]
		if resampling == 0:
			s0*=scale
		else:
			s0=resample(s0,resampling,chord_i=chord_i,reversed=isReversed[i])*scale

		if transl != None:
			s0=s0+transl

		s1=extractPointsFromSTEPFile(path+inFiles[i+1])#[1:-1,:]


		if type(mirrorCoord) == int:
			s1[:,mirrorCoord]*=-1.

		if type(mirrorCoord) == list:
			for m in mirrorCoord:
				s1[:,m]*=-1.

		if reverse[i+1]:
			s1=s1[::-1]
		if resampling == 0:
			s1*=scale
		else:
			s1=resample(s1,resampling,chord_i=chord_i,reversed=isReversed[i])*scale

		if transl != None:
			s1+=transl

		geometry.addSection(s0,s1)
		geometry.sections[len(geometry.sections)-1].setNumberOfIterpolants(nIterps[i])
	nSec=len(geometry.getInterpolants())/2
	return geometry,nSec

def makePanels(pr1,wing,wakeVector,nWakePanels,wingGeometry,nSec,wakeGeometry=None,swapWake=True):
	geometrySection={}
	wingInterpolants=wingGeometry.getInterpolants()
	wing.nSpanwise=len(wingInterpolants)-1

	if wakeGeometry is not None:
		wakeInterpolants=wakeGeometry.getInterpolants()

	for i in range(0,wing.nSpanwise):
		oldPanels=[wing.panels[k] for k in wing.panels]

		rSec=wingInterpolants[i]
		lSec=wingInterpolants[i+1]
		if wakeGeometry == None:
			makeSpanwiseTwoLoftFoilReg(pr1.dom1,wing,rSec,lSec,wakeVector=wakeVector,nWakePanels=nWakePanels,swapWake=swapWake)
		else:
			rwSec=wakeInterpolants[i]
			lwSec=wakeInterpolants[i+1]
			makeSpanwiseTwoLoftFoilReg(pr1.dom1,wing,rSec,lSec,wakeVector=None,nWakePanels=0,wakeA=rwSec,wakeB=lwSec,swapWake=swapWake)

		geometrySection[i]=PanelGroup()
		newPanels=[wing.panels[k] for k in wing.panels]
		addedPanels=[item for item in newPanels if item not in oldPanels]
		for j in range(0,len(addedPanels)):
			geometrySection[i].addPanel(addedPanels[j])

	nPan=len(addedPanels)
	for i in range(0,wing.nSpanwise):
		s=geometrySection[i]
		for p in s.getPanelList():
			pan=s.getPanel(p)
			if i >0 :
				pan.connectivity['d4']=pan.gid-nPan
			if i <(wing.nSpanwise-1):
				pan.connectivity['d2']=pan.gid+nPan
	return geometrySection

def extractLDM(REG,alpha):
	f0=REG.getForceAtRefPT()
	f=f0
	f=f0
	cl=-f[0]*np.sin(alpha)+f[1]*np.cos(alpha)
	cd=f[0]*np.cos(alpha)+f[1]*np.sin(alpha)
	cm=f[5]
	return cl,cd,cm

def reportRegion(fname,r):
	f=open(fname,'w')
	for i in range(0,r.nPanels):
		p=r.panels[i]
		f.write(str(i)+'\t'+str(p.gid)+'\t'+str(p.pressureCoef)+'\t'+str(p.G)+'\t'+str(p.S)+'\t'+str(np.linalg.norm(p.forceVector))+'\n')
	f.close()
