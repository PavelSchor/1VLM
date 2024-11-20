import numpy as np

def proj(u,v):
	return u*np.dot(v,u)/np.dot(u,u)

def grammSmithOrthoNormalize(V):
	U=np.copy(V)
	for i in range(1,V.shape[0]):
		for j in range(0,i):
			U[i] -= proj(U[j],V[i])
		U[i] /= np.sqrt( (U[i]**2).sum(axis=0))
	return U

def _getTransformMatrix(targetCSYS,sourceCSYS):
	cs1=np.zeros(sourceCSYS.shape)
	cs2=np.zeros(sourceCSYS.shape)
	res=np.zeros(sourceCSYS.shape)
	for i in range(0,sourceCSYS.shape[0]):
		cs2[i]=sourceCSYS[i]/np.linalg.norm(sourceCSYS[i])
		cs1[i]=targetCSYS[i]/np.linalg.norm(targetCSYS[i])
	for i in range(0,cs2.shape[0]):
		for j in range(0,cs2.shape[1]):
			res[i][j]=np.dot(cs1[i],cs2[j])

	return res

def getDCMFromCYSES(targetCSYS,sourceCSYS):
	return _getTransformMatrix(targetCSYS,sourceCSYS)


#void vect_getTransformMatrix(double*res,double*targetCSYS,double*sourceCSYS){
#	// directional cosine matrix DCM
#	int i,j;
#	double vm1,vm2,cs1[9],cs2[9];
#
#	for(i=0;i<3;i++){
#		vm1=vect_magn(&targetCSYS[i*3]);
#		vm2=vect_magn(&sourceCSYS[i*3]);
#		for(j=0;j<3;j++){
#			cs1[i*3+j]=targetCSYS[i*3+j]/vm1;
#			cs2[i*3+j]=sourceCSYS[i*3+j]/vm2;
#		}
#	}
#	for(i=0;i<3;i++){
#		for(j=0;j<3;j++){
#			res[i*3+j]= vect_dot(&cs1[i*3],&cs2[j*3]);







def transformVector(vect,targetCSYS,sourceCSYS):
	return np.dot( _getTransformMatrix(targetCSYS,sourceCSYS),vect)

def makeCSYS_from2Vectors(u,v,ax=1):
	ax=int(ax)
	n=u.shape[0]
	res=np.zeros((n,n))
	res[0]=u/np.linalg.norm(u)
	if ax == 1:
		res[1]=v/np.linalg.norm(v)
		res[2]=np.cross(res[0],res[1])
	if ax ==2:
		res[2]=v/np.linalg.norm(v)
		res[1]=-np.cross(res[0],res[2])
	return grammSmithOrthoNormalize(res)

def getDCM(inp):
	Phi,Theta,Psi=inp[0],inp[1],inp[2]
	#	Euler Angles:
	#		Phi,	Roll Angle, rad
	#		Theta,	Pitch Angle, rad
	#		Psi,	Yaw Angle, rad

	sinR = np.sin(Phi);
	cosR = np.cos(Phi);
	sinP = np.sin(Theta);
	cosP = np.cos(Theta);
	sinY = np.sin(Psi);
	cosY = np.cos(Psi);

	#H(1,1) = cosP * cosY;
	#H(1,2) = cosP * sinY;
	#H(1,3) = -sinP;
	#H(2,1) = sinR * sinP * cosY - cosR * sinY;
	#H(2,2) = sinR * sinP * sinY + cosR * cosY;
	#H(2,3) = sinR * cosP;
	#H(3,1) = cosR * sinP * cosY + sinR * sinY;
	#H(3,2) = cosR * sinP * sinY - sinR * cosY;
	#H(3,3) = cosR * cosP;

	return 	np.array([ 	[cosP * cosY, cosP * sinY,  -sinP],
				[sinR * sinP * cosY - cosR * sinY,  sinR * sinP * sinY + cosR * cosY, sinR * cosP],
				[ cosR * sinP * cosY + sinR * sinY,  cosR * sinP * sinY - sinR * cosY, cosR * cosP]
				])

def rotateAxisAngle(inp,ax,r,org=None):
	#http://answers.google.com/answers/threadview/id/361441.html
	#http://www.cprogramming.com/tutorial/3d/quaternions.html
	if org is None:
		org=np.zeros(3)
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
		
def rotateCSYSByAngles(csys,angles):
	csn=np.copy(csys)
	for i in range(0,3):
		csn=rotateAxisAngle(csn,csys[i],angles[i])
		
	return csn