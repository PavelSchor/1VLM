from utils.mesh import *
from vlm.vlmMethod1 import *
import numpy as np
np.setbufsize(8192*2)
import pylab as plt
import time
from tqdm import tqdm




cg=np.array([0.25,0,0])
chord=1.
span=10.
alpha=np.deg2rad(10)
beta=np.deg2rad(0)
sp=35.
q=0.5*1*sp**2
VEL=np.array([np.cos(alpha),-np.sin(beta),np.sin(alpha)])*sp
pr1=VlmProblem()
self=pr1.dom1
pr1.dom1.COR=cg

xDir=np.array([1,0,0])
yDir=np.array([0,1,0])
zDir=np.array([0,0,1])

wing_l=PanelGroup()
nc=21
pts1=np.zeros((nc,3)); pts1[:,0]=np.linspace(0,1,nc);pts1[:,1]=-5.
pts2=np.zeros((nc,3)); pts2[:,0]=np.linspace(0,1,nc);pts2[:,1]=0.
pts3=np.zeros((nc,3)); pts3[:,0]=np.linspace(0,1,nc);pts3[:,1]=5.
ll=[pts1,pts2,pts3]
nip=[10,10]
makeMeshFromArrayList(pr1,wing_l,ll,nip=nip)
appendLoftedWakeRegKnownLasts(pr1,wing_l,1,wakeVector=np.array([30,0,0]))

pr1.addRegion(wing_l)
pr1.dom1.rho=1.
pr1.dom1.initSolution()
pr1.dom1.VEL=VEL


alphaRng=np.deg2rad(np.linspace(-10,20,31))


CL=np.zeros_like(alphaRng)
CD=np.zeros_like(alphaRng)
CM=np.zeros_like(alphaRng)

start = time.time()
for i in range(0,len(alphaRng)):
	alpha=alphaRng[i]
	VEL=np.array([np.cos(alpha),0.1,np.sin(alpha)])*sp
	pr1.dom1.VEL=VEL


	#%time pr1.dom1.compute()
	#%time pr1.dom1.computeViscous()
	#%timeit pr1.dom1.compute()


	##%timeit 
	#pr1.dom1.computeViscous(10)
	
	#pr1.dom1.compute()
	#%time 
	#pr1.dom1.updateTimeStep()
	pr1.dom1.compute()
	#print abs(pr1.dom1.RESIDUAL).max()
	CF=self.FORCE.sum(axis=0)/q/chord/span
	CL[i]=-CF[0]*np.sin(alpha)+CF[2]*np.cos(alpha)
	CD[i]= CF[0]*np.cos(alpha)+CF[2]*np.sin(alpha)
	CM[i]=self.MOMENT[1]/q/chord/span/chord
	#pr1.plotToFile('figs/wing_plain{:d}.vtp'.format(i),True)
	
plt.plot(np.rad2deg(alphaRng),CL);plt.show();
plt.plot(np.rad2deg(alphaRng),CM);plt.show();
