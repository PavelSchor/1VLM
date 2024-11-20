import numpy as np
import stl
#from stl import mesh

def extractPtsFromSTL(fname):
	m=stl.mesh.Mesh.from_file(fname)
	return np.vstack((m.v0,m.v1,m.v2))
