#%%

from src import data_postprocess


#%%

data = np.c_[ecm.time, ecm.voltage, vt]
data_reduced = data_postprocess.data_desample(data, x=0, y=1, stepsize=0.001)

f = open("expval/Exp_gitt_terminvolt.txt", "w")
# np.savetxt(f, np.c_[ecm.time[0:-1:20], ecm.voltage[0:-1:20], vt[0:-1:20]], fmt='%1.6e', delimiter=', ')
np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("expval/Exp_gitt_OCVpts.txt", "w")
np.savetxt(f, np.c_[ecm.time[ecm.restEnd], ecm.voltage[ecm.restEnd]], fmt='%1.6e', delimiter=', ')
f.close()


#%%
f = open("expval/Exp_gitt_rmse.txt", "w")
no = np.arange (1,26)
np.savetxt(f, np.c_[no, rmse_ecm], fmt='%12.4f', delimiter=', ')
f.close()

#%%
t_curve = t_curve - t_curve[0]
data = np.c_[t_curve[0:2300], v_curve[0:2300], vt[0:2300]]
data_reduced = data_postprocess.data_desample(data, x=0, y=1, stepsize=0.015)
#%%
f = open("expval/layered_method_expvoltb.txt", "w")
# np.savetxt(f, np.c_[t_curve[0:-1:10], v_curve[0:-1:10], vt[0:-1:10]], fmt='%1.6e', delimiter=', ')
np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()


#%%

data = np.c_[ecm.time, zx_all, soc]
data_reduced = data_postprocess.data_desample(data, x=0, y=10, stepsize=0.002)

#%
f = open("expval/Exp_gitt_soc.txt", "w")
# np.savetxt(f, np.c_[ecm.time[0:-1:50], zx_all[0:-1:50,:], soc[0:-1:50]], fmt='%1.6e', delimiter=', ')
np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()



#%%

data = np.c_[t_curve[0:3690:1], zx_all[0:3690:1], z_curve[0:3690:1]]
data_reduced = data_postprocess.data_desample(data, x=0, y=10, stepsize=0.01)

f = open("expval/Exp_gitt_soc_5th.txt", "w")
# np.savetxt(f, np.c_[t_curve[0:3690:10], zx_all[0:3690:10], z_curve[0:3690:10]], fmt='%1.6e', delimiter=', ')
np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()

#%%
# data = np.c_[t_curve[0:2082:1], Ix_all[0:2082:1], i_curve[0:2082:1]]
# data_reduced = data_postprocess.data_desample(data, x=0, y=10, stepsize=0.01)

f = open("expval/Exp_gitt_I_5th.txt", "w")
np.savetxt(f, np.c_[t_curve[0:2082:10], Ix_all[0:2082:10], i_curve[0:2082:10]], fmt='%1.6e', delimiter=', ')
# np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()


#%%
f = open("expval/Exp_gitt_timescale_NE.txt", "w")
np.savetxt(f, np.c_[soc_mid, timescale_5, timescale_10, timescale_15, timescale_20], fmt='%1.6e', delimiter=', ')
f.close()


#%%

err = np.abs(ecm.voltage - vt)
data = np.c_[soc, ecm.voltage, vt, err]
data_trim = data_postprocess.data_trim(data)
data_reduced = data_postprocess.data_desample(data_trim)

#%% 
# 0.4C CC predication
f = open("expval/Exp_CC_0.4C.txt", "w")
np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()
#%%
# 2C CC predication
f = open("expval/Exp_CC_2C.txt", "w")
err = np.abs(ecm.voltage - vt)
np.savetxt(f, data_reduced, fmt='%1.6e', delimiter=', ')
f.close()

#%%
# 1C CC predication
f = open("expval/Exp_CC_0.5C.txt", "w")
err = np.abs(ecm.voltage - vt)
inter = 40
np.savetxt(f, np.c_[soc[0:-1:inter], ecm.voltage[0:-1:inter], vt[0:-1:inter], err[0:-1:inter]], fmt='%1.6e', delimiter=', ')
f.close()


#%% 
# parameters

data =  np.c_[soc_mid, params_array[:, [2,3]]]
f = open("expval/Exp_gitt_params.txt", "w")
np.savetxt(f, data, fmt='%1.6e', delimiter=', ')
f.close()
#%%
data =  np.c_[soc_plot, Rd, R0]
f = open("expval/Exp_gitt_params_polyfit.txt", "w")
np.savetxt(f, data, fmt='%1.6e', delimiter=', ')
f.close()






