import numpy as np
from scipy.interpolate import griddata
import pylab as plt
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import interp2d

class AeroSurfFactory(object):
	def __init__(self):
		self.alphaRng=np.linspace(-10,20,51)
		self.pts_a=[]
		self.pts_y=[]
		self.pts_r=[]
		self.xi=np.zeros((2,2))
		
	def setAlphaRange(self,arng):
		self.alphaRng=arng

	def setYi(self,yi):
		self.yi=yi
	

	def addAirfoilData(self,ypos,a1,interpolatingLiftCurve=None,order=1):
		#if interpolatingLiftCurve != None:
##			af0=np.interp(interpolatingLiftCurve[:,1],a1[:,0],a1[:,1])
##			af1=np.interp(self.alphaRng,interpolatingLiftCurve[:,0],af0)
			#spl = UnivariateSpline(interpolatingLiftCurve[:,1],interpolatingLiftCurve[:,0])
			#ai=spl(a1[:,0])

			#s = InterpolatedUnivariateSpline(ai,a1[:,1], k=order)
			#af1=s(self.alphaRng)

		#else:
##			af1=np.interp(self.alphaRng,a1[:,0],a1[:,1])
			#s = InterpolatedUnivariateSpline(a1[:,0],a1[:,1], k=order)
			#af1=s(self.alphaRng)
		af1=ll=np.interp(self.alphaRng,a1[:,0],a1[:,1])

		idx=np.abs(self.yi-ypos).argmin()
		y=self.yi[idx]*np.ones(len(self.alphaRng))
		self.pts_a.extend(self.alphaRng.tolist())
		self.pts_y.extend(y.tolist())
		self.pts_r.extend(af1.tolist())

	def make(self):
		X,Y=np.meshgrid(self.yi,self.alphaRng)
		self.pts_i=np.array([self.pts_y,self.pts_a]).T
		vals=np.array(self.pts_r)
		Z=griddata(self.pts_i,vals,(X,Y),method='linear')
		self.X,self.Y,self.Z=X,Y,Z
		#self.interpolator = interp2d(self.X, self.Y, self.Z, kind='linear')
		
	def getData(self,xi):#yi,alpha):
		#return self.interpolator(yi,alpha)
		#return griddata(self.pts_i,self.pts_r,np.hstack((yi[None].T,alpha[None].T)),method='linear')
		return griddata(self.pts_i,self.pts_r,xi,method='linear')
		
	def plot(self):
		fig = plt.figure()
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(self.X,self.Y,self.Z,rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
		plt.show()

	def save(self,fname):
		np.savetxt(fname,self.Z.flatten(order='f'))

class WingFactory(object):
	def __init__(self):
		self.n=31
		self.span=10.
		self.chordsTable=np.array([ [-5e32,1],[5e32,1] ])
		self.twistsTable=np.array([ [-5e32,0],[5e32,0] ])
		self.liftData=AeroSurfFactory()
		self.dragData=AeroSurfFactory()
		self.momentData=AeroSurfFactory()
		self.dir_main='data/'
		self.dir_geom=self.dir_main+'geometry/'
		self.dir_lift=self.dir_main+'lift/'

	def setN(self,n):
		if ( n % 2 == 0.):
			self.n=n+1
		else:
			self.n=n

	def setSpan(self,b):
		self.span=b/2.0
		self.fi=np.zeros(self.n)
		for i in range(0,self.n):
			self.fi[i]=np.pi*(i+1)/(self.n+1);
		self.yi=-self.span*np.cos(self.fi)
		self.liftData.setYi(self.yi)
		self.dragData.setYi(self.yi)
		self.momentData.setYi(self.yi)

	def getYi(self):
		self.setSpan(self.span)
		return self.yi

	def makeAeroData(self):
		self.liftData.make()
		self.dragData.make()
		self.momentData.make()

	def save(self):
		self.makeAeroData()
		np.savetxt(self.dir_geom+'yi.txt',self.yi*-1.)
		np.savetxt(self.dir_geom+'chords.txt',np.interp(self.yi,self.chordsTable[:,0],self.chordsTable[:,1]))
		np.savetxt(self.dir_geom+'wing_twist.txt',np.interp(self.yi,self.twistsTable[:,0],self.twistsTable[:,1]))
		np.savetxt(self.dir_geom+'span.txt',np.array([self.span]))
		np.savetxt(self.dir_lift+'alpha.txt',np.deg2rad(self.liftData.alphaRng))
		self.liftData.save(self.dir_lift+'lift_surf.txt')
		self.dragData.save(self.dir_lift+'drag_surf.txt')
		self.momentData.save(self.dir_lift+'moment_surf.txt')		
		


