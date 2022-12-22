import numpy as np

import sys
from src import shell_vol_ratios
from scipy import optimize
from scipy import interpolate


class shellECM:
	"""
	"""
	def __init__(self, data, allids, params):
		"""
		Initialize with HPPC battery cell data and model parameters.
		"""
		self.time = data.time
		# convert current sign to positive for discharging
		self.current = data.current * -1.0
		self.voltage = data.voltage

		self.initiPt = allids[0]
		self.discSta = allids[1]
		self.discEnd = allids[2]
		self.restSta = allids[3]
		self.restEnd = allids[4]

		self.Q_cell = params.Q_cell

	def soc(self):
		"""
		"""
		current = self.current
		time = self.time

		Q_cell = self.Q_cell
		dt = np.diff(time)

		nc = len(current)
		z = np.ones(nc)

		for k in range(1, nc):

			# current step
			# z[k] = z[k-1] - current[k] * dt[k-1] / Q_cell

			# previous step
			# z[k] = z[k-1] - current[k-1] * dt[k-1] / Q_cell

			# previous + current step
			z[k] = z[k - 1] - (current[k] +
							   current[k - 1]) / 2 * dt[k - 1] / Q_cell

		return z

	def ocv(self, soc, pts=False, vz_pts=None):
		"""
		"""
		if pts is True:

			id4 = self.restEnd
			v_pts = np.append(self.voltage[0], self.voltage[id4])
			z_pts = np.append(soc[0], soc[id4])
			# i_pts = np.append(self.current[0], self.current[id4])
			# t_pts = np.append(self.time[0], self.time[id4])

			# -----------------------------------------------------------
			# extrapolate to 0 so that [0,1]
			z_diffs = abs(np.diff(z_pts))
			npts_extra = np.ceil(
				(min(z_pts) - 0.0) / (np.mean(z_diffs))).astype('int')
			# print(npts_extra)
			z_pts_extra = np.linspace(min(z_pts), 0.0, npts_extra + 1)

			slope_1 = (v_pts[-2] - v_pts[-1]) / (z_pts[-2] - z_pts[-1])
			slope_2 = (v_pts[-3] - v_pts[-2]) / (z_pts[-3] - z_pts[-2])
			slope_3 = (v_pts[-4] - v_pts[-3]) / (z_pts[-4] - z_pts[-3])
			delta_slope = 0.5 * (slope_1 - slope_2) + 0.5 * (slope_2 - slope_3)

			# -----------------------------------------------------------
			# way 1
			# slope = slope_1 + delta_slope
			# v_pts_extra = v_pts[-1] - (z_pts[-1] - z_pts_extra) * slope

			# or way 2
			# v_pts_extra = np.ones(npts_extra+1) * v_pts[-1]
			# for i in range(1, npts_extra+1):
			#     v_pts_extra[i] = v_pts_extra[i-1] - (z_pts_extra[i-1] - z_pts_extra[i]) * (slope_1 + i * delta_slope)

			# or way 3
			# adjust them against pseudo OCV
			# x = np.array([z_pts[-1], z_pts[-2], z_pts[-3]])
			# y = np.array([v_pts[-1], v_pts[-2], v_pts[-3]])
			x = z_pts[-1:-4:-1]
			y = v_pts[-1:-4:-1]
			p = np.poly1d(np.polyfit(x, y, 2))
			v_pts_extra = p(z_pts_extra)
			# -----------------------------------------------------------

			# remove the first that duplicate min(z_pts) / z_pts[-1]
			z_pts = np.append(z_pts, z_pts_extra[1:])
			v_pts = np.append(v_pts, v_pts_extra[1:])

			# -----------------------------------------------------------
			# add more pointd below 0 and above 1
			z_pts_below = -0.01
			v_pts_below = p(z_pts_below)

			z_pts_above = 1.05
			slope = (v_pts[1] - v_pts[0]) / (z_pts[1] - z_pts[0])
			v_pts_above = v_pts[0] + slope * (z_pts_above - z_pts[0])

			z_pts = np.append(z_pts, z_pts_below)
			v_pts = np.append(v_pts, v_pts_below)

			z_pts = np.insert(z_pts, 0, z_pts_above)
			v_pts = np.insert(v_pts, 0, v_pts_above)

			#-----------------------------------------------------------
			# linear interpolation
			#-----------------------------------------------------------
			# # f = interpolate.interp1d(z_pts, v_pts)
			# # ocv = f(soc)

			ocv = np.interp(soc, z_pts[::-1], v_pts[::-1])

			#-----------------------------------------------------------
			# polyfit interpolation
			#-----------------------------------------------------------
			# p = np.poly1d( np.polyfit(z_pts, v_pts, 8) )
			# ocv = p(soc)

			# return ocv, i_pts, t_pts, v_pts, z_pts
			return ocv, v_pts, z_pts

		elif vz_pts is not None:

			# v_pts, z_pts = vz_pts
			v_pts = vz_pts[:,0]
			z_pts = vz_pts[:,1]

			# # f = interpolate.interp1d(z_pts, v_pts)
			# # ocv = f(soc)
			ocv = np.interp(soc, z_pts[::-1], v_pts[::-1])

			#-----------------------------------------------------------
			# polyfit interpolation
			#-----------------------------------------------------------
			# p = np.poly1d( np.polyfit(z_pts, v_pts, 8) )
			# ocv = p(soc)

			return ocv
		else:
			print('wrong input parameters!')

	@staticmethod
	def extend_array(soc, rx_array):
		"""
		extend r0_array to the same length of soc
		"""
		#-------------------------------------------------------
		# method 1
		z_pts = np.ravel(rx_array[:, [0, 1]])
		r_pts = np.repeat(rx_array[:, 2], 2)

		rx_array_ext = np.interp(soc, z_pts[::-1], r_pts[::-1])

		#-------------------------------------------------------
		# method 2
		# z_pts = (rx_array[:,0] + rx_array[:,1]) / 2
		# r_pts = rx_array[:, 2]

		# z_pts = np.insert(z_pts, 0, 1.0)
		# z_pts = np.append(z_pts, 0.0)

		# r_pts = np.insert(r_pts, 0, rx_array[0, 2])
		# r_pts = np.append(r_pts, rx_array[-1, 2])

		# rx_array_ext = np.interp(soc, z_pts[::-1], r_pts[::-1])

		# f = interpolate.interp1d(z_pts, r_pts)
		# rx_array_ext = f(soc)

		#-------------------------------------------------------
		# method 3

		# z_pts = (rx_array[:,0] + rx_array[:,1]) / 2
		# r_pts = rx_array[:, 2]

		# p = np.poly1d( np.polyfit(z_pts, r_pts, 2) )
		# rx_array_ext = p(soc)

		#-------------------------------------------------------

		return rx_array_ext


	def get_R0(self, soc):

		id0 = self.initiPt
		id1 = self.discSta
		id2 = self.discEnd
		id3 = self.restSta
		id4 = self.restEnd

		nrow = len(self.initiPt)

		R0_array = np.zeros((nrow, 3))

		for k in range(0, 25):

			start = id0[k]
			end = id4[k]

			# calculate R_0
			di = abs(self.current[id1[k]] - self.current[id0[k]])
			dv = abs(self.voltage[id1[k]] - self.voltage[id0[k]])
			r0_down = dv / di
			di = abs(self.current[id3[k]] - self.current[id2[k]])
			dv = abs(self.voltage[id3[k]] - self.voltage[id2[k]])
			r0_up = dv / di
			r0 = np.mean(np.array([r0_down, r0_up]))

			R0_array[k, 0] = soc[start]
			R0_array[k, 1] = soc[end]
			R0_array[k, 2] = r0_down

		return R0_array



	def check_some_pulse(self, soc, vz_pts, para_array, nr_layer, i_pulse):

		r_diff = para_array[0]
		id0 = self.initiPt
		# id1 = self.discSta
		# id2 = self.discEnd
		id3 = self.restSta
		id4 = self.restEnd

		k = i_pulse - 1
		start = id0[k]
		end = id4[k]
		t_curve = self.time[start:end + 1]
		i_curve = self.current[start:end + 1]
		v_curve = self.voltage[start:end + 1]
		z_curve = soc[start:end + 1]

		staPt = id0[k] - start
		endPt = id4[k] - start

		zx_all, Ix_all = self.multi_volt_sources_pulse(t_curve, i_curve, z_curve,
												 vz_pts, r_diff, nr_layer)

		z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]

		ocp = self.ocv(z_surf, False, vz_pts)

		vt = ocp - para_array[1] * i_curve

		return t_curve, i_curve, v_curve, z_curve, ocp, vt, zx_all, Ix_all

	def get_paras_shellECM(self, soc, vz_pts, paras_ini, nr_layer):
		"""
		non_least_squares
		https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html

		"""

		nrow = len(self.initiPt)

		# first two columns includes the soc range
		# soc_sta, soc_end, r_diff, r0
		paras_array = np.zeros((nrow, 4))

		id0 = self.initiPt
		id1 = self.discSta
		id2 = self.discEnd
		id3 = self.restSta
		id4 = self.restEnd

		for k in range(0, 25):

			start = id0[k]
			end = id4[k]
			t_curve = self.time[start:end + 1]
			i_curve = self.current[start:end + 1]
			v_curve = self.voltage[start:end + 1]
			z_curve = soc[start:end + 1]

			# calculate R_0
			di = abs(self.current[id1[k]] - self.current[id0[k]])
			dv = abs(self.voltage[id1[k]] - self.voltage[id0[k]])
			r0_down = dv / di
			di = abs(self.current[id3[k]] - self.current[id2[k]])
			dv = abs(self.voltage[id3[k]] - self.voltage[id2[k]])
			r0_up = dv / di
			r0 = np.mean(np.array([r0_down, r0_up]))

			staPt = id0[k] - start
			endPt = id4[k] - start

			bound_range = ([0.0, 0.0], [10, 10])
			paras_ini[1] = r0

			result = optimize.least_squares(self.fun_obj_nls,
											paras_ini,
											bounds=bound_range,
											args=(t_curve, i_curve, z_curve,
												  v_curve, vz_pts, staPt,
												  endPt, nr_layer),
											verbose=1)

			paras_array[k, 0] = soc[start]
			paras_array[k, 1] = soc[end]
			paras_array[k, 2:] = result.x

		return paras_array

	def fun_obj_nls(self, fit_paras, t_curve, i_curve, z_curve, v_curve,
					vz_pts, staPt, endPt, nr_layer):

		r_diff = fit_paras[0]

		zx_all, Ix_all = self.multi_volt_sources_pulse(t_curve, i_curve, z_curve,
												 vz_pts, r_diff, nr_layer)

		z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]

		ocp = self.ocv(z_surf, False, vz_pts)

		vt = ocp - fit_paras[1] * i_curve

		return vt[staPt:endPt + 1] - v_curve[staPt:endPt + 1]


	def multi_volt_sources_pulse(self, time, current, soc, vz_pts, r_diff_1,
						   nr_layer):

		# note current flowing out a particle as positive
		# exp data discharging as negative while positive in pybamm

		rat_R, rat_V = shell_vol_ratios.get_ratios(nr_layer)

		R_diffs = np.ones(nr_layer - 1)
		# r_diff_1---the innermost layer diffusion resistance
		for i in range(0, nr_layer - 1):
			R_diffs[i] = r_diff_1 * rat_R[i]

		if not all(R_diffs > 0.0):
			sys.exit('internal resistance not ALL positive!')

		Q_cell = self.Q_cell
		dt = np.diff(time)
		time_steps = len(time)

		# initiate concentration
		# [c_1, c_2, c_3, ..., c_n, c_N]
		z_all = np.zeros((time_steps, nr_layer))
		z_all[0, :] = np.ones((1, nr_layer)) * soc[0]

		# initiate current
		# [I1, I2, ..., In, I_N]
		I_all = np.zeros((time_steps, nr_layer))

		# loop over time starting from the first time step
		# to calculate variables at the second time step
		# the first timestep hold initial values
		for k in range(0, time_steps - 1):

			# the driving force is concentration gradient at last timestep
			Volt_all = z_all[k, :]

			# KCL to get I bar, gradient driven mass flux
			I_bar = np.zeros(nr_layer)
			for i in range(0, nr_layer - 1):
				I_bar[i] = (Volt_all[i] - Volt_all[i + 1]) / R_diffs[i]
			I_bar[nr_layer - 1] = current[k]

			# update current in shell brunches
			# particle center bc: null flux
			I_all[k, 0] = I_bar[0].copy()
			for i in range(1, nr_layer):
				I_all[k, i] = I_bar[i] - I_bar[i - 1]

			z_all[k + 1, :] = z_all[k, :] - I_all[k, :] * dt[k] / (Q_cell * rat_V)

		return z_all, I_all


	def multi_volt_sources(self, soc, r_diff_array, nr_layer):

		rat_R, rat_V = shell_vol_ratios.get_ratios(nr_layer)

		r_diff_array_ext = self.extend_array(soc, r_diff_array)

		current = self.current
		time = self.time
		Q_cell = self.Q_cell
		dt = np.diff(time)
		time_steps = len(time)

		# initiate concentration
		# [c_1, c_2, c_3, ..., c_n, c_N]
		z_all = np.zeros((time_steps, nr_layer))
		z_all[0, :] = np.ones((1, nr_layer)) * soc[0]

		# initiate current
		# [I1, I2, ..., In, I_N]
		I_all = np.zeros((time_steps, nr_layer))

		R_diffs = np.ones(nr_layer - 1)

		for k in range(0, time_steps - 1):

			for i in range(0, nr_layer - 1):
				R_diffs[i] = r_diff_array_ext[k] * rat_R[i]
				# R_diffs[i] = r_diff_array * rat_R[i]

			# the driving force is concentration gradient at last timestep
			Volt_all = z_all[k, :]

			# KCL to get I bar, gradient driven mass flux
			I_bar = np.zeros(nr_layer)
			for i in range(0, nr_layer - 1):
				I_bar[i] = (Volt_all[i] - Volt_all[i + 1]) / R_diffs[i]
			I_bar[nr_layer - 1] = current[k]

			# update current in shell brunches
			# particle center bc: null flux
			I_all[k, 0] = I_bar[0].copy()
			for i in range(1, nr_layer):
				I_all[k, i] = I_bar[i] - I_bar[i - 1]

			z_all[k + 1, :] = z_all[k, :] - I_all[k, :] * dt[k] / (Q_cell * rat_V)

		return z_all, I_all



	def get_vt(self, soc, vz_pts, params_array, nr_layer):
		"""
		Determine voltage from equivalent circuit model.
		"""

		r_diff_array = params_array[:,[0, 1, 2]]
		zx_all, Ix_all = self.multi_volt_sources(soc, r_diff_array, nr_layer)

		z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]
		ocp = self.ocv(z_surf, False, vz_pts)

		r0_array = params_array[:, [0, 1, 3]]
		r0_array_ext = self.extend_array(soc, r0_array)

		vt = ocp - r0_array_ext * self.current

		return vt


# -------------------------------------------------------------------------------
	def get_paras_shellECM_diff(self, soc, vz_pts, diff_ini, nr_layer):
		"""
		non_least_squares
		https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html

		"""

		nrow = len(self.initiPt)

		# first two columns includes the soc range
		# soc_sta, soc_end, r_diff, r0
		paras_array = np.zeros((nrow, 4))

		id0 = self.initiPt
		id1 = self.discSta
		id2 = self.discEnd
		id3 = self.restSta
		id4 = self.restEnd

		for k in range(0, 25):

			start = id0[k]
			end = id4[k]
			t_curve = self.time[start:end + 1]
			i_curve = self.current[start:end + 1]
			v_curve = self.voltage[start:end + 1]
			z_curve = soc[start:end + 1]

			# calculate R_0
			di = abs(self.current[id1[k]] - self.current[id0[k]])
			dv = abs(self.voltage[id1[k]] - self.voltage[id0[k]])
			r0_down = dv / di
			di = abs(self.current[id3[k]] - self.current[id2[k]])
			dv = abs(self.voltage[id3[k]] - self.voltage[id2[k]])
			r0_up = dv / di
			r0 = np.mean(np.array([r0_down, r0_up]))
			# r0 = r0_up

			staPt = id2[k] - start
			endPt = id4[k] - start

			bound_range = ([0.0, 10.0])

			result = optimize.least_squares(self.fun_obj_nls_diff,
											diff_ini,
											bounds=bound_range,
											args=(t_curve, i_curve, z_curve,
												  v_curve, vz_pts, staPt,
												  endPt, r0, nr_layer),
											verbose=1)

			paras_array[k, 0] = soc[start]
			paras_array[k, 1] = soc[end]
			paras_array[k, 2] = result.x
			paras_array[k, 3] = r0

		return paras_array

	def fun_obj_nls_diff(self, fit_paras, t_curve, i_curve, z_curve, v_curve,
					vz_pts, staPt, endPt, r0, nr_layer):

		r_diff = fit_paras

		zx_all, Ix_all = self.multi_volt_sources_pulse(t_curve, i_curve, z_curve,
												 vz_pts, r_diff, nr_layer)

		z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]

		ocp = self.ocv(z_surf, False, vz_pts)

		vt = ocp - r0 * i_curve

		return vt[staPt:endPt + 1] - v_curve[staPt:endPt + 1]


# ---------------------------------------------------------------------------------




# -------------------------------------------------------------------------------
	def get_paras_shellECM_uni(self, soc, vz_pts, paras_ini, nr_layer):
		"""
		non_least_squares
		https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html

		"""

		id0 = self.initiPt
		id1 = self.discSta
		id2 = self.discEnd
		id3 = self.restSta
		id4 = self.restEnd

		start = id0[0]
		end = id2[0]
		bound_range = ([0.0, 0.0], [10, 10])

		result = optimize.least_squares(self.fun_obj_nls_uni,
										paras_ini,
										bounds=bound_range,
										args=(soc, vz_pts, nr_layer, start, end),
										verbose=1)

		# paras_array = result.x

		return result.x

	def fun_obj_nls_uni(self, fit_paras, soc, vz_pts, nr_layer, start, end):

		r_diff = fit_paras[0]

		zx_all, Ix_all = self.multi_volt_sources_uni(soc, r_diff, nr_layer)

		z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]

		ocp = self.ocv(z_surf, False, vz_pts)

		vt = ocp - fit_paras[1] * self.current

		return vt[start:end+1] - self.voltage[start:end+1]


	def multi_volt_sources_uni(self, soc, r_diff, nr_layer):

		rat_R, rat_V = shell_vol_ratios.get_ratios(nr_layer)

		current = self.current
		time = self.time
		Q_cell = self.Q_cell
		dt = np.diff(time)
		time_steps = len(time)

		# initiate concentration
		# [c_1, c_2, c_3, ..., c_n, c_N]
		z_all = np.zeros((time_steps, nr_layer))
		z_all[0, :] = np.ones((1, nr_layer)) * soc[0]

		# initiate current
		# [I1, I2, ..., In, I_N]
		I_all = np.zeros((time_steps, nr_layer))

		R_diffs = np.ones(nr_layer - 1)

		for k in range(0, time_steps - 1):

			for i in range(0, nr_layer - 1):
				R_diffs[i] = r_diff * rat_R[i]

			# the driving force is concentration gradient at last timestep
			Volt_all = z_all[k, :]

			# KCL to get I bar, gradient driven mass flux
			I_bar = np.zeros(nr_layer)
			for i in range(0, nr_layer - 1):
				I_bar[i] = (Volt_all[i] - Volt_all[i + 1]) / R_diffs[i]
			I_bar[nr_layer - 1] = current[k]

			# update current in shell brunches
			# particle center bc: null flux
			I_all[k, 0] = I_bar[0].copy()
			for i in range(1, nr_layer):
				I_all[k, i] = I_bar[i] - I_bar[i - 1]

			z_all[k + 1, :] = z_all[k, :] - I_all[k, :] * dt[k] / (Q_cell * rat_V)

		return z_all, I_all

# -------------------------------------------------------------------------------