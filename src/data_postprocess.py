
import numpy as np

def data_trim(data):

	n_pts = len(data[:,0])

	for k in range(1, n_pts):

		if data[k,1] < 2.5:
			data_trim = data[0:k,:]
			break

	return data_trim




def data_desample(data, x = 0, y = 1, stepsize = 0.01):

	# step size to control data density 
	# stepsize = 0.01

	# data = np.c_[soc, ecm.voltage]

	# number of data points
	n_pts = len(data[:,x])

	x_scale = np.max(data[:,x]) - np.min(data[:,x])
	y_scale = np.max(data[:,y]) - np.min(data[:,y])

	data_reduced = []
	data_reduced.append(data[0,:])

	ref_pts = 0 

	stop = 0 
	for i in range(0, int(1e8)):

		for k in range(ref_pts+1, n_pts-1):

			dx = (data[k,x] - data[ref_pts,x]) / x_scale
			dy = (data[k,y] - data[ref_pts,y]) / y_scale
			dis = np.sqrt(dx ** 2 + dy ** 2)

			if dis > stepsize:
				data_reduced.append(data[k,:])
				ref_pts = k
				break

			if k == n_pts - 2:
				data_reduced.append(data[n_pts - 1,:])
				ref_pts = k + 1
				stop = 1
				break

		if stop:
			break

	return np.array(data_reduced)