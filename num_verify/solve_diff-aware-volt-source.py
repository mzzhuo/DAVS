
# ----------------------------------------------
# solve the diffusion-aware voltage source 
# produce results in Section 3.1
# ----------------------------------------------

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

# from src import params_plett
from src import params_chen2020PE
# from src import shell_vol_ratios
from src import shellVoltSource

#%%
# read file
# filepath = "test/j_p_6s.txt"
filepath = "test/current_6s.txt"
data = pd.read_csv(filepath)

#%%
nr_layer = 20

# shellvs = shellVoltSource(data, params, nr_layer)
sign = -1.0
shellvs = shellVoltSource(data, params_chen2020PE, nr_layer, sign)

#%%

time_f, c_all_f, I_all_f = shellvs.solve_forward_Euler()

#%%
time_b, c_all_b, I_all_b = shellvs.solve_backward_Euler()


#%% surface concentration and ocp

c_p_s = 1.5*c_all_b[:,-1]-0.5*c_all_b[:,-2]
sto = c_p_s / params_chen2020PE.c_max
ocp = ( -0.8090 * sto + 4.4875 
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )


fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(time_b, ocp, 'r-', label='Shell-imp')
ax1.legend()
plt.show()

f = open("test/OCP_p_6s-ecm.txt", "w")
np.savetxt(f, np.c_[time_b, ocp], fmt='%1.6e', delimiter=', ')
f.close()


#%% comaper with plett
filepath = "test/c_s_6s.txt"
data_pb = pd.read_csv(filepath)
time_pb = data_pb['time']
cs_pb = data_pb['cs']

#%%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
# ax1.plot(time_pb, cs_pb/1000, 'k-', label='FVM')
ax1.plot(time_f, (1.5*c_all_f[:,-1]-0.5*c_all_f[:,-2])/1000, 'b-', label='Shell-exp')
ax1.plot(time_b, (1.5*c_all_b[:,-1]-0.5*c_all_b[:,-2])/1000, 'r-', label='Shell-imp')
ax1.legend()
plt.show()

# %%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(time_b, c_all_b[:,0]/1000, 'b-', label='Shell-exp')
ax1.plot(time_b, c_all_b[:,10]/1000, 'r-', label='Shell-imp')
ax1.plot(time_b, c_all_b[:,-1]/1000, 'r-', label='Shell-imp')
ax1.legend()
plt.show()

#%%
r_p = np.linspace(0.025, 0.975, 20)

f = open("test/c_p_6s-ecm.txt", "w")
# np.savetxt(f, np.c_[time_b, c_all_b], fmt='%1.6e', delimiter=', ')
np.savetxt(f, np.c_[r_p, np.transpose(c_all_b)], fmt='%1.6e', delimiter=', ')

f.close()

#%%
f = open("test/c_p_s_6s-ecm.txt", "w")
np.savetxt(f, np.c_[time_b, (1.5*c_all_b[:,-1]-0.5*c_all_b[:,-2])], fmt='%1.6e', delimiter=', ')
f.close()




# %%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(c_all_b[:,0]/1000, 'r-', label='c_surf')
ax1.plot(c_all_b[:,10]/1000, 'g-', label='c_surf')
ax1.plot(c_all_b[:,-1]/1000, 'b-', label='c_surf')
ax1.legend()
plt.show()

# %%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
# ax1.plot(c_all_b[:,0]/1000, 'r-', label='c_surf')
# ax1.plot(c_all_b[:,10]/1000, 'g-', label='c_surf')
ax1.plot(c_all_b[100,:]/1000, 'b-', label='c_surf')
ax1.legend()
plt.show()

#%% comaper with plett
filepath = "test/cse_plett.txt"
data = np.loadtxt(filepath)
time = np.arange(0,len(data)) / 3600

# %%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(time, c_all_f[:,-1]/1000, 'r-', label='Shell-exp')
ax1.plot(time[0:-1], c_all_b[:,-1]/1000, 'b-', label='Shell-imp')
ax1.plot(time, data/1000, 'k-', label='FVM')
ax1.legend()
plt.xlabel("time (h)", fontsize=16)
plt.ylabel("concentration (kmol.m-3)", fontsize=16)
plt.show()