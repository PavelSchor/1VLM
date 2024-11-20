#
def polyArea(poly):
	A=0.0
	for i in range(0,poly.shape[0]-1):
		A+= (poly[i,0]*poly[i+1,1] - poly[i+1,0]*poly[i,1])
		return 0.5*A

