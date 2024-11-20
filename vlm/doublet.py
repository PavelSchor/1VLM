import numpy as np

def VXJE_doublet_edge_influence(node_a, node_b, x):
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

	print 'z,m,e1,e2,r1,r2   ' , z,m,e1,e2,r1,r2
	print 'u,v  ' ,u,v
	if (u == v):
		delta_theta = 0.0;
	else:
		delta_theta = np.arctan2(u - v, 1 + u * v);

	res = delta_theta;

	return res;

