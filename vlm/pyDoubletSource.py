import numpy as np

def doubletEdgeInfl(x, node_a, node_b):
	d = np.sqrt(pow(node_b[0] - node_a[0], 2) + pow(node_b[1] - node_a[1], 2));
	if np.abs(d) < 1e-10:
		return 0.0
	z = x[2];
	m = (node_b[1] - node_a[1]) / (node_b[0] - node_a[0]);
	e1 = pow(x[0] - node_a[0], 2) + pow(z, 2);
	e2 = pow(x[0] - node_b[0], 2) + pow(z, 2);
	r1 = np.sqrt(e1 + pow(x[1] - node_a[1], 2));
	r2 = np.sqrt(e2 + pow(x[1] - node_b[1], 2));
	h1 = (x[0] - node_a[0]) * (x[1] - node_a[1]);
	h2 = (x[0] - node_b[0]) * (x[1] - node_b[1]);
	u = (m * e1 - h1) / (z * r1);
	v = (m * e2 - h2) / (z * r2);
	if np.all(u == v):
		delta_theta = 0.0;
	else:
		delta_theta = np.arctan2(u - v, 1 + u * v);
	print "doublet edge: ",node_a,node_b,delta_theta
	return delta_theta

def sourceEdgeInfl(x, node_a, node_b):
	d = np.sqrt(pow(node_b[0] - node_a[0], 2) + pow(node_b[1] - node_a[1], 2));
	if np.abs(d) < 1e-10:
		return 0.0
	z = x[2];
	m = (node_b[1] - node_a[1]) / (node_b[0] - node_a[0]);
	e1 = pow(x[0] - node_a[0], 2) + pow(z, 2);
	e2 = pow(x[0] - node_b[0], 2) + pow(z, 2);
	r1 = np.sqrt(e1 + pow(x[1] - node_a[1], 2));
	r2 = np.sqrt(e2 + pow(x[1] - node_b[1], 2));
	h1 = (x[0] - node_a[0]) * (x[1] - node_a[1]);
	h2 = (x[0] - node_b[0]) * (x[1] - node_b[1]);
	u = (m * e1 - h1) / (z * r1);
	v = (m * e2 - h2) / (z * r2);
	if np.all(u == v):
		delta_theta = 0.0;
	else:
		delta_theta = np.arctan2(u - v, 1 + u * v);
        return ((x[0] - node_a[0]) * (node_b[1] - node_a[1]) - (x[1] - node_a[1]) * (node_b[0] - node_a[0])) / d * np.log((r1 + r2 + d) / (r1 + r2 - d)) - np.abs(z) * delta_theta; 

class pyDSPanel(object):
	def __init__(self):
		self.G=0.0
		self.S=0.0
		pass

	def setPoints(self,pts):
		self.pts=pts
		self.ptsl=self.transform(self.orig,self.pts)

	def setCSYS(self,csys):
		self.csys=csys
	
	def setOrig(self,orig):
		self.orig=orig

	def transform(self,orig,pts):
#		ptsl=np.zeros(pts.shape)
		ptsl=pts-orig
		if len(ptsl.shape) >1:
			for i in range(0,ptsl.shape[0]):
				ptsl[i]=np.dot(self.csys,ptsl[i])
		if len(ptsl.shape) == 1:
			ptsl=np.dot(self.csys,ptsl)
		return ptsl

	def setG(self,G):
		self.G=G

	def setS(self,S):
		self.S=S
	
	def getDoubletInfluence(self,PT):
		PTL=self.transform(self.orig,PT)
		infl=0.0
		for i in range(0,self.ptsl.shape[0]):
			infl+=doubletEdgeInfl(PTL,self.ptsl[i] , self.ptsl[(i+1)%self.ptsl.shape[0]])
		return infl

	def getSourceInfluence(self,PT):
		PTL=self.transform(self.orig,PT)
		infl=0.0
		for i in range(0,self.ptsl.shape[0]):
			infl+=sourceEdgeInfl(PTL,self.ptsl[i] , self.ptsl[(i+1)%self.ptsl.shape[0]])
		return infl

	def getInfluence(self,PT,G=None,S=None):
		if G==None:
			G=self.G
		if S==None:
			S=self.S
		infl=0.0
		infl+= G/4./np.pi*(self.getDoubletInfluence(PT))
		infl+= S/4./np.pi*(self.getSourceInfluence(PT))
		return infl
