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
from utils.bCad import *
#from am_analysis.amFunctions import *


#from fem.fem3 import FEM1
#from fem.beam3d import Beam3D
from utils.interactor import PanelBeamInteractor
from utils.loftedSurface import *
from utils.elasticWing import ElasticWing1
import utils.csys as cs
from utils.beamCalc import BeamCrossSection as BeamCalc
from fem.utils.writeBDF import BDFWriter
from utils.simpleSTEPReader import extractPointsFromSTEPFile
from utils.geometryGen import *
from flightLoads.pointLoads import *

import affinity
import multiprocessing
 
affinity.set_process_affinity_mask(0,2**multiprocessing.cpu_count()-1)

bdfWriter=BDFWriter()
vect=B3vect()
bc=BCad()

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def semiCosineSpacing(nt,nc,rc):
	nr=nt-nc
	r=np.zeros(nt)
	r[0:nc]=(1+np.cos(np.linspace(np.pi,np.pi/2,nc)))*rc
	r[nc-1:nt]=np.logspace(np.log10(r[nc-1]), np.log10(1.0), num=nr+1)
	
	#r[nc-1:nt]=np.linspace(rc,1,nr+1)
	return r

def semiCosineSpacingString(nt,nc,rc,scl=1):
	r=semiCosineSpacing(nt,nc,rc)*scl
	s=''
	for i in range(0,len(r)):
		s+=str(r[i])+' '
	return s

def exponentialSpacing(nt,nc=1e3,base=10):
	nr=nt-nc
	r=np.zeros(nt)
	r=(np.logspace(np.log10(1), np.log10(nc+1), num=nt)-1)/nc
	#r[0:nc]=(np.logspace(np.log10(1), np.log10(1001), num=nc,base=base)-1)*1e-3*rc
	return r


def panuklChordwiseSpacingString(nt,nc,rc,scl=1,mode='exponential'):
	if mode=='semiCosine':
		r=semiCosineSpacing(nt,nc,rc)*scl
	if mode=='exponential':
		r=exponentialSpacing(nt,nc)*scl
	s=''
	for i in range(0,len(r)):
		s+=str(r[i])+' '
	return s


def extractSubsectionFromPanuklFile(fname,gid0=0,makeSymmetricSections=[]):
	panels={}
	#f=open('pokus2.inp')
	f=open(fname)
	ll=f.readlines()
	f.close()
	pg={}
	lc=0
	try:
		npan=int(ll[lc].split()[0])
		nwpan1=int(ll[lc].split()[1])
		nwpan=int(ll[lc].split()[2])
	except:
		lc=1
		npan=int(ll[lc].split()[0])
		nwpan1=int(ll[lc].split()[1])
		nwpan=int(ll[lc].split()[2])
		
	def panuklPGGridToPanels(ll,lc,panels,gid0):
		n=int(ll[lc].split()[0])
		for i in range(0,n):
			
			lc+=1
			l=ll[lc]
			#print l
			if len(l.split()) == 1:
				break
			x=np.array(l.split()[1::],dtype=np.float).reshape((4,3))
			panels[len(panels)]=doubletSourcePanel(x[1],x[0],x[3],x[2])
			panels[len(panels)-1].gid=int(l.split()[0])-1+gid0
		return lc+1

	def isConnectiblePanelID(x):
		if x==0 or x==100000:
			return False
		else:
			return True
		


	#lc=0
	lc=panuklPGGridToPanels(ll,lc,panels,gid0)
	#while len(ll[lc].split()) == 1:
		#lc=panuklPGGridToPanels(ll,lc,panels)
	lc+=nwpan+1
	
	for i in range(0,len(panels)):
		l=ll[lc]
		if len(l.split())<2:
			break
		ii=int(l.split()[0])-1
		
		if ii == 0:
			print ii, l
		x=np.array(l.split()[1::],dtype=np.int)
		#x-=1
		p=panels[ii]
		
		try:
			if isConnectiblePanelID(x[1]):
				p.connectivity['d3']=x[1]-1+gid0
		except:
				pass
		try:
			if isConnectiblePanelID(x[5]):
				p.connectivity['d1']=x[5]-1+gid0
		except:
				pass
		try:
			if isConnectiblePanelID(x[3]):
				p.connectivity['d2']=x[3]-1+gid0
		except:
				pass
		try:
			if isConnectiblePanelID(x[7]):
				p.connectivity['d4']=x[7]-1+gid0
		except:
				pass
		for conn in p.connectivity:
			try:
				p.connectivity[conn]=abs(p.connectivity[conn])
			except:
				pass
		panels[i].gid=int(l.split()[0])-1+gid0
		if ii == 0:
			print ii, l
			print panels[i].gid,panels[i].connectivity

		lc+=1

	lc+=nwpan1
	#print ll[lc]
	nPG=int(ll[lc])
	for i in range(0,nPG):
		lc+=1
		l=ll[lc]
		start=int(l.split()[0])
		stop=int(l.split()[1])+1
		name=str(l.split()[2])
		if not name in makeSymmetricSections:
			pg[name]=subSection()
			for j in range(start,stop):
				pg[name].addPanel(panels[j])
		else:
			pg[name+'_1']=subSection()
			pg[name+'_2']=subSection()
			for j in range(start,start+(stop-start)/2):
				pg[name+'_1'].addPanel(panels[j])
			for j in range(start+(stop-start)/2,stop):
				pg[name+'_2'].addPanel(panels[j])

	return pg

def extractSubsectionFromPanukl2017File(fname,gid0=0):
	panels={}
	#f=open('pokus2.inp')
	f=open(fname)
	ll=f.readlines()
	f.close()
	pg={}
	lc=0
	if 'ver' in ll[0].split():
		lc=1
		
	try:
		npan=int(ll[lc].split()[0])
		nwpan1=int(ll[lc].split()[1])
		nwpan=int(ll[lc].split()[2])
	except:
		lc=1
		npan=int(ll[lc].split()[0])
		nwpan1=int(ll[lc].split()[1])
		nwpan=int(ll[lc].split()[2])
		
	def panuklPGGridToPanels(ll,lc,panels,gid0):
		n=int(ll[lc].split()[0])
		for i in range(0,n):
			
			lc+=1
			l=ll[lc]
			#print l
			if len(l.split()) == 1:
				break
			x=np.array(l.split()[1::],dtype=np.float).reshape((4,3))
			panels[len(panels)]=doubletSourcePanel(x[1],x[0],x[3],x[2])
			panels[len(panels)-1].gid=int(l.split()[0])-1+gid0
		return lc+1

	def isConnectiblePanelID(x):
		if x==0 or x==100000:
			return False
		else:
			return True
		


	#lc=0
	lc=panuklPGGridToPanels(ll,lc,panels,gid0)
	#while len(ll[lc].split()) == 1:
		#lc=panuklPGGridToPanels(ll,lc,panels)
	lc+=nwpan+1
	
	for i in range(0,len(panels)):
		l=ll[lc]
		if len(l.split())<2:
			break
		ii=int(l.split()[0])-1
		
		#if ii == 0:
			#print ii, l
		x=np.array(l.split()[1::],dtype=np.int)
		x=np.abs(x)
		p=panels[ii]
		p.gid=int(l.split()[0])-1+gid0
		try:
			if isConnectiblePanelID(x[1]):
				p.connectivity['d3']=x[1]-1+gid0
		except:
				pass
		try:
			if isConnectiblePanelID(x[5]):
				p.connectivity['d1']=x[5]-1+gid0
		except:
				pass
		try:
			if isConnectiblePanelID(x[3]):
				p.connectivity['d2']=x[3]-1+gid0
		except:
				pass
		try:
			if isConnectiblePanelID(x[7]):
				p.connectivity['d4']=x[7]-1+gid0
		except:
				pass
		#for conn in p.connectivity:
			#try:
				#p.connectivity[conn]=abs(p.connectivity[conn])
			#except:
				#pass
		
		#if ii == 0:
			#print ii, l
			#print panels[i].gid,panels[i].connectivity

		lc+=1

	lc+=nwpan1
	#print ll[lc]
	nPG=int(ll[lc])
	print nPG
	for i in range(0,nPG):
		lc+=1
		l=ll[lc]
		start=int(l.split()[0])
		stop=int(l.split()[1])+1
		name=str(l.split()[9])
		pg[name]=subSection()
		print start, stop, name
		for j in range(start,stop):
			pg[name].addPanel(panels[j])
	return pg



def makePanuklWingToFile1(fname,p1,p2,airfoils,ns,incidence={'ax':np.array((0,0,0)),'orig':np.array((0,0,0)),'r':0.0},dihedral={'ax':np.array((0,0,0)),'orig':np.array((0,0,0)),'r':0.0}):
	le=np.copy(p1)
	dl=np.linalg.norm(p2-p1,axis=1)
	pp1=np.copy(p1)
	pp2=np.copy(p2)
	if not incidence['r'] == 0.0:
		pp1=bc.rotAx(p1,incidence['ax'],incidence['orig'],incidence['r'])
		pp2=bc.rotAx(p2,incidence['ax'],incidence['orig'],incidence['r'])
		print incidence['r']
	if not dihedral['r'] == 0.0:
		pp1=bc.rotAx(pp1,dihedral['ax'],dihedral['orig'],dihedral['r'])
		pp2=bc.rotAx(pp2,dihedral['ax'],dihedral['orig'],dihedral['r'])
		print dihedral['r']
	dp=pp2-pp1
	twist=-np.arctan2(dp[:,2],dp[:,0])

	sc=0
	f=open(fname,'w')
	for i in range(0,len(airfoils)):
		sc+=ns[i]

		f.write(airfoils[i]+'\n')
		f.write('{:6.6f}\n'.format(dl[i]))
		f.write('{:6.6f}\t{:6.6f}\t{:6.6f}\n'.format(le[i,0],abs(le[i,1]),le[i,2]))
		f.write('{:6.6f}\t{:6.6f}\t{:6.6f}\n'.format(np.rad2deg(dihedral['r']),np.rad2deg(twist[i]),0.0))
		f.write('{:d}\n'.format(int(sc)))
	
	f.close()
		#/home/pavel/PanuklProjects/profile/n64212a1.prf
#0.89785
#5.56900	0.88000	1.55389
#0.00000	-3.50000	0.00000
#10


def makeSparPointsLeTePct(le,te,pct):
	v=te-le
	dl=np.linalg.norm(v,axis=1)
	return le+v/dl*pct

def makePanuklPrfToFile1(fname,name,fnt,fnb,ix=0,iy=2,length='abs'):
	tp=extractPointsFromSTEPFile(fnt)
	bp=extractPointsFromSTEPFile(fnb)
	l0=100.
	if length=='abs':
		scl=l0/np.linalg.norm(tp[-1]-tp[0])
	if length=='x':
		scl=l0/(tp[-1][ix] - tp[0][ix])
	tp-=tp[0]
	bp-=bp[0]
	tp*=scl
	bp*=scl
	
	n=len(tp)
	
	f=open(fname,'w')
	f.write(' '+str(n)+'\t '+name+'\n')
	for i in range(0,n-1):
		f.write('{:8.8f} {:8.8f} {:8.8f} {:8.8f}\n'.format(tp[i,ix],tp[i,iy],bp[i,ix],bp[i,iy]))
	i+=1
	f.write('{:8.8f} {:8.8f} {:8.8f} {:8.8f}\n'.format(l0,tp[i,iy],l0,bp[i,iy]))
	f.close()
