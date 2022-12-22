# ----------------------------------------------
# parameterize the shell ECM in Section 4.1 
# generate prediction results of the shell ECM in Section 4.2
# ----------------------------------------------


import matplotlib.pyplot as plt

import numpy as np
# import pandas as pd
# import codecs

from src import pulseData
# from src import shellVoltSource
from src import shellECM
from src import pseudoOCV
from src import params_chen2020Cell
from src import params_chen2020PE
from src import params_chen2020NE


# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# 
# GITT 4% pulse parameterization
# 
# ----------------------------------------------------------------------------------------------------


#%%
file = './NDK01-67 - CC, 4pc SOC pulse, 20pc SOC pulse - 1C - cell 2.mpt'

#%%
data_1 = pulseData.process_one_cycle(file, 1)

# find the 5 ids of each pulse
# allids = pulseData.get_indices_discharge(data_1)
allids = data_1.get_indices_discharge()

# remove redundant data before the first pulse anneDrive/shellECN/icell/pulse_shell.py')
# update all the ids
allids = data_1.trim(allids)

#%%
# fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
# ax1.plot(data_1.time, data_1.voltage, 'ro', markersize=1, label='voltage')
# ax1.plot(data_1.time[allids[0]], data_1.voltage[allids[0]], 'o', markersize=4, label='voltage')
# ax1.plot(data_1.time[allids[1]], data_1.voltage[allids[1]], 'x', markersize=4, label='voltage')
# ax1.plot(data_1.time[allids[2]], data_1.voltage[allids[2]], '>', markersize=4, label='voltage')
# ax1.plot(data_1.time[allids[3]], data_1.voltage[allids[3]], '<', markersize=4, label='voltage')
# ax1.plot(data_1.time[allids[4]], data_1.voltage[allids[4]], '+', markersize=4, label='voltage')


#%%
# from src import shell_vol_ratios
# rat_R, rat_V = shell_vol_ratios.get_ratios(10)

# from src import diffResRatios
# rat_R2, rat_q, _ = diffResRatios.get_drRatio(10, 'eq_thickness')

#%%

ecm = shellECM(data_1, allids, params_chen2020Cell)

soc = ecm.soc()


ocv, v_pts, z_pts = ecm.ocv(soc, pts=True)
vz_pts=np.c_[v_pts, z_pts]

# soc_sm = np.linspace(1.0, 0.0, 100)
# ocv_sm = ecm.ocv(soc_sm, False, vz_pts)

# np.save('expval/vzpts.npy', vz_pts)

#%%
# ocv_vol, ocv_soc = pseudoOCV.get_ocv()
# vz_pts=np.c_[ocv_vol, ocv_soc]


#%% 
# # output cell soc and volt
# SoC_cell = np.linspace(1.0, 0.0, 400)
# OCV_cell = ecm.ocv(SoC_cell, False, vz_pts)
# # Q_cell = SoC_cell * 5

# f = open("expval\OCV_cell2.csv", "w")
# np.savetxt(f, np.c_[SoC_cell, OCV_cell], fmt='%1.6e', delimiter=', ')
# f.close()

#%%
fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
ax1.plot(soc, ocv, 'b-', linewidth=1, markersize=2, label='ocv')
ax1.plot(z_pts, v_pts, 'ro', markersize=4, label='ocv_pts')
# ax1.plot(ocv_soc, ocv_vol, 'r', markersize=4, label='ocv_pts')

# ax1.plot(soc_sm, ocv_sm, 'k-', linewidth=1, markersize=2, label='ocv')
# ax1.plot(SoC_cell, OCV_cell, 'k', markersize=4, label='ocv_pts')

ax1.legend(fontsize=14)
ax1.invert_xaxis()

plt.show()






#%%
# nr_layer = 10
# paras_ini = np.array([1.0, 30e-3])
# params_array = ecm.get_paras_shellECM(soc, vz_pts, paras_ini, nr_layer)


#%%
nr_layer = 10

diff_ini = 1.0
params_array = ecm.get_paras_shellECM_diff(soc, vz_pts, diff_ini, nr_layer)


#%%
nr_layer = 10

# np.save('expval/params_' + str(nr_layer) + '_layers', params_array)
params_array = np.load('expval/params_' + str(nr_layer) + '_layers-sep.npy')
# vz_pts = np.load('expval/vzpts.npy')

#%%
# get R0 array
# R0_array = ecm.get_R0(soc)

#%%
soc_mid = (params_array[:,0] + params_array[:,1]) / 2
timescale = 4 * np.pi * params_chen2020PE.radius ** 3 * 96485.3321 * params_chen2020PE.c_max / 10 * params_array[:,2] * params_chen2020PE.N

fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
ax1.plot(soc_mid, timescale/3600, 'o-', markersize=4)
plt.show()


#%%

# nr_layer = 5
# params_array_5 = np.load('expval/params_' + str(nr_layer) + '_layers-sep.npy')

# nr_layer = 10
# params_array_10 = np.load('expval/params_' + str(nr_layer) + '_layers-sep.npy')

# nr_layer = 15
# params_array_15 = np.load('expval/params_' + str(nr_layer) + '_layers-sep.npy')

# nr_layer = 20
# params_array_20 = np.load('expval/params_' + str(nr_layer) + '_layers-sep.npy')

# soc_mid = (params_array_10[:,0] + params_array_10[:,1]) / 2

# # coeff = 4 * np.pi * params_chen2020PE.radius ** 3 * 96485.3321 * params_chen2020PE.c_max * params_chen2020PE.N / 3600
# coeff = 4 * np.pi * params_chen2020NE.radius ** 3 * 96485.3321 * params_chen2020NE.c_max * params_chen2020NE.N / 3600

# timescale_5 = coeff / 5 * params_array_5[:,2]
# timescale_10 = coeff / 10 * params_array_10[:,2]
# timescale_15 = coeff / 15 * params_array_15[:,2]
# timescale_20 = coeff / 20 * params_array_20[:,2]

# fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
# ax1.plot(soc_mid, timescale_5, 'o-', markersize=4)
# ax1.plot(soc_mid, timescale_10, 's-', markersize=4)
# ax1.plot(soc_mid, timescale_15, 's-', markersize=4)
# ax1.plot(soc_mid, timescale_20, 'o-', markersize=4)
# plt.show()

# coeff_PE = 4 * np.pi * params_chen2020PE.radius ** 3 * params_chen2020PE.c_max * params_chen2020PE.N 
# coeff_NE = 4 * np.pi * params_chen2020NE.radius ** 3 * params_chen2020NE.c_max * params_chen2020NE.N 

#%%




#%%
soc_mid = (params_array[:,0] + params_array[:,1]) / 2
# fitting 
p = np.poly1d( np.polyfit(soc_mid, params_array[:,3], 2) )
soc_plot = np.linspace(0.0, 1.0, 100)
R0 = p(soc_plot)

fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
ax1.plot(soc_mid, params_array[:,3], 'o', markersize=4)
ax1.plot(soc_plot, R0, '-', markersize=4)
plt.show()

#%%
soc_mid = (params_array[:,0] + params_array[:,1]) / 2
# fitting 
p = np.poly1d( np.polyfit(soc_mid, params_array[:,2], 2) )
soc_plot = np.linspace(0.0, 1.0, 100)
Rd = p(soc_plot)

fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
ax1.plot(soc_mid, params_array[:,2], 'o', markersize=4)
ax1.plot(soc_plot, Rd, '-', markersize=4)
plt.show()






#%%
# nr_layer = 10
r_diff_array = params_array[:,[0, 1, 2]]
# r_diff_array = np.mean( params_array[:,2])
zx_all, Ix_all = ecm.multi_volt_sources(soc, r_diff_array, nr_layer)

#%%
# vt = ecm.get_vt(soc, vz_pts, params_array, nr_layer)
z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]
ocp = ecm.ocv(z_surf, False, vz_pts)
r0_array = params_array[:, [0, 1, 3]]
r0_array_ext = ecm.extend_array(soc, r0_array)
# r0_array_ext = np.mean( params_array[:,3])

vt = ocp - r0_array_ext * ecm.current
        

#%%
fig, ax1 = plt.subplots(tight_layout=True)

ax1.plot(ecm.time, ecm.voltage, 'b-', markersize=1, label='exp')
ax1.plot(ecm.time, vt, 'r-', markersize=1, label='ecm')
# ax1.plot(ecm.time, np.abs(ecm.voltage-vt), 'b-', markersize=1, label='ocp')

ax1.legend()
plt.show()


# %%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(ecm.time, ecm.current, 'k-', label='curr')

for i in range(0,nr_layer):
    La = 'I' + str(i+1)
    ax1.plot(ecm.time, Ix_all[:,i], '-', lw=1, label=La)

ax1.legend()
plt.show()



# %%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)

for i in range(0,nr_layer):
    La = 'soc' + str(i+1)
    ax1.plot(ecm.time, zx_all[:,i], '-', lw=1, label=La)

ax1.plot(ecm.time, soc, 'k-', lw=1.5, label='bulk soc')
# ax1.plot(ecm.time[0::100], z2[0::100], 'C5-', label='z2')

ax1.set_xlabel(r'Time $t$ (s)',fontsize=12 )
ax1.set_ylabel(r'SoC',fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=10) 

ax1.legend()
plt.show()



#%%
# calculate rmse
nrow = len(ecm.initiPt)

id0 = ecm.initiPt
id3 = ecm.restSta
id4 = ecm.restEnd

rmse_ecm = np.zeros((nrow, 1))

for k in range(0, nrow):

    start = id3[k]
    end = id4[k]

    vt_pulse_exp = ecm.voltage[start:end+1]
    vt_pulse_ecm = vt[start:end+1]

    rmse_ecm[k] = np.linalg.norm( vt_pulse_ecm - vt_pulse_exp ) / np.sqrt( len(vt_pulse_exp) )
#%
fig, ax1 = plt.subplots(figsize=(8, 5), tight_layout=True)

ax1.plot(rmse_ecm, 'C4s--', markersize=3, lw=1, label='shellECM')

ax1.set_xlabel(r'pulse number',fontsize=12 )
ax1.set_ylabel(r'rmse (V)',fontsize=12)
ax1.tick_params(axis='both', which='major', labelsize=10)

ax1.legend(fontsize=14)

plt.show()



#%%
nr_layer = 10
ipse = 14
# para_array_p = np.array([50e-2, 25e-3])
paras_array_p = params_array[ipse-1,2:]
#%
# t_curve, i_curve, v_curve, z_curve, v1_pulse, vt_pulse, zx_all, Ix_all = ecm.check_some_pulse(soc, vz_pts, para_array, 10, ipse)
t_curve, i_curve, v_curve, z_curve, ocp, vt, zx_all, Ix_all = ecm.check_some_pulse(soc, vz_pts, paras_array_p, nr_layer, ipse)

#%%
fig, ax1 = plt.subplots(tight_layout=True)
# t_curve = t_curve - t_curve[0]

ax1.plot(t_curve, z_curve, 'k-', markersize=1, label='exp')
for i in range(0,nr_layer):
    La = 'z' + str(i+1)
    ax1.plot(t_curve, zx_all[:,i], '-', label=La)

ax1.legend()
plt.show()

#%%
fig, ax1 = plt.subplots(tight_layout=True)

ax1.plot(t_curve, i_curve, 'k-', markersize=1, label='exp')
for i in range(0,nr_layer):
    La = 'I' + str(i+1)
    ax1.plot(t_curve, Ix_all[:,i], '-', label=La)

ax1.legend()
plt.show()



#%%
fig, ax1 = plt.subplots(tight_layout=True)
# t_curve = t_curve - t_curve[0]

ax1.plot(t_curve, ocp, 'r-', markersize=1, label='exp')
ax1.plot(t_curve, v_curve, 'k-', markersize=1, label='exp')
ax1.plot(t_curve, vt, 'b-', markersize=1, label='exp')

ax1.legend()
plt.show()


#%%
lay_cen = np.linspace(0.05, 0.95, 10)

t_curve = t_curve - t_curve[0]
#%%
time_plot = [0, 30, 60, 100, 150, 200, 300, 700, 2500]
idx = []
for itime in time_plot:
    for k in range(1, len(t_curve)):
        if abs(t_curve[k] - itime) < 0.3:
            idx.append(k)
            break
#%
fig, ax1 = plt.subplots(figsize=(8, 5), tight_layout=True)
for i in range(0,len(idx)):
    La = str( ecm.time[ idx[i] ].astype('int') ) + ' s'
    ax1.plot(lay_cen, zx_all[idx[i],:], '-o', markersize=4, label=La)

ax1.set_xlim(0, 1)
# ax1.set_ylim(0.0,1.05)

ax1.legend(fontsize=14)

ax1.legend()
ax1.set_xlabel(r'layer center',fontsize=14 )
ax1.set_ylabel(r'soc',fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=12)

plt.show()








# ----------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------
# 
# CC predication
# 
# ----------------------------------------------------------------------------------------------------

#%%

# file = 'D:/Clouds/OneDrive - Imperial College London/shellECN/Nialls data/CC discharge + pulse data/1C and 2C data/2C/NDK01-67 - CC, 4pc SOC pulse, 20pc SOC pulse - 2C - cell 1.mpt'
file = './NDK01-30 - b2 - CC, OCV+P(4pc), OCV+P(20pc) at  0,4C_02_MB_CB2.mpt'



#%%
data_0 = pulseData.process_one_cycle(file, 0)

# find the 5 ids of each pulse
# allids = pulseData.get_indices_discharge(data_1)
allids = data_0.get_indices_discharge()

# remove redundant data before the first pulse and 
# update all the ids
allids = data_0.trim(allids)

# #%%
# fig, ax1 = plt.subplots(figsize=(8, 6), tight_layout=True)
# ax1.plot(data_0.time, data_0.voltage, 'ro', markersize=1, label='voltage')
# ax1.plot(data_0.time[allids[0]], data_0.voltage[allids[0]], 'o', markersize=4, label='voltage')
# ax1.plot(data_0.time[allids[1]], data_0.voltage[allids[1]], 'x', markersize=4, label='voltage')
# ax1.plot(data_0.time[allids[2]], data_0.voltage[allids[2]], '>', markersize=4, label='voltage')
# ax1.plot(data_0.time[allids[3]], data_0.voltage[allids[3]], '<', markersize=4, label='voltage')
# ax1.plot(data_0.time[allids[4]], data_0.voltage[allids[4]], '+', markersize=4, label='voltage')

#%%
ecm = shellECM(data_0, allids, params_chen2020Cell)
soc = ecm.soc()

#%%
r_diff_array = params_array[:,[0, 1, 2]]
zx_all, Ix_all = ecm.multi_volt_sources(soc, r_diff_array, nr_layer)

#%%
# vt = ecm.get_vt(soc, vz_pts, params_array, nr_layer)
z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]
ocp = ecm.ocv(z_surf, False, vz_pts)

#%%
r0_array = params_array[:, [0, 1, 3]]
r0_array_ext = ecm.extend_array(soc, r0_array)
# r0_array_ext = np.mean( params_array[:,3])
vt = ocp - r0_array_ext * ecm.current


#%%
# ocv_vol, ocv_soc = pseudoOCV.get_ocv()
# vz_pts=(ocv_vol, ocv_soc)
# #%%
# nr_layer = 10

# paras_ini = np.array([1.0, 25e-3])
# params_array = ecm.get_paras_shellECM_uni(soc, vz_pts, paras_ini, nr_layer)

# #%%
# r_diff = params_array[0]
# zx_all, Ix_all = ecm.multi_volt_sources_uni(soc, r_diff, nr_layer)

# z_surf = 1.5 * zx_all[:, -1] - 0.5 * zx_all[:, -2]
# ocp = ecm.ocv(z_surf, False, vz_pts)

# r0 = params_array[1]
# vt = ocp - r0 * ecm.current




#%%
fig, ax1 = plt.subplots(tight_layout=True)

# ax1.plot(ecm.time, vt, 'r-', markersize=1, label='ecm')
# ax1.plot(ecm.time, ecm.voltage, 'bo', markersize=1, label='exp')
ax1.plot(ecm.time, vt-ecm.voltage, 'bo', markersize=1, label='exp')

# ax1.plot(ecm.time, ocp, 'k-', markersize=1, label='ocp')

ax1.legend()
plt.show()
