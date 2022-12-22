# ----------------------------------------------
# solve the shell ECM containing diffusion-aware voltage source 
# reproduce SPM simulation results
# produce results in Section 3.3
# ----------------------------------------------


import numpy as np
import pandas as pd

from src import params_chen2020PE
from src import params_chen2020NE
from src import shellVoltSource

from scipy import constants



#%%
# read applied current data
# filepath = "test/current.txt"
# filepath = "test/current_gitt.txt"
# filepath = "test/current_CC_halfC.txt"
filepath = "test/current_CC_3C.txt"

data = pd.read_csv(filepath)

#%%
nr_layer = 20

# applied current positive for discharging
# both PE NE particle lithium leaving as positive
sign = -1.0
shellvs_pe = shellVoltSource(data, params_chen2020PE, nr_layer, sign)


sign = 1.0
shellvs_ne = shellVoltSource(data, params_chen2020NE, nr_layer, sign)


#%%
time, c_all_ne, I_all_ne = shellvs_ne.solve_backward_Euler()

#%%
time, c_all_pe, I_all_pe = shellvs_pe.solve_backward_Euler()


#%% 
R_cons = constants.R
F_cons = constants.physical_constants["Faraday constant"][0]
T_cons = 298.15
V_ref = R_cons * T_cons / F_cons

#%%
# PE
# calculate particle surface concentration
c_p_surf = 1.5 * c_all_pe[:,-1] - 0.5 * c_all_pe[:,-2]
# calculate exchange current density
j0_p = 3.42e-6 * 1000 ** 0.5 * c_p_surf ** 0.5 * (params_chen2020PE.c_max - c_p_surf) ** 0.5
# calculate the charge transfer resistance
# linear approximation
# R_ct_p = V_ref / j0_p / (params_chen2020PE.V_ed * params_chen2020PE.epsilon * params_chen2020PE.Sa_vol)  
# related to applied current
I = 15
R_ct_p = 2 * V_ref / I * np.arcsinh(I / (params_chen2020PE.V_ed * params_chen2020PE.epsilon * params_chen2020PE.Sa_vol) / 2 / j0_p);

#%%
# NE
c_n_surf = 1.5 * c_all_ne[:,-1] - 0.5 * c_all_ne[:,-2]
j0_n = 6.48e-7 * 1000 ** 0.5 * c_n_surf ** 0.5 * (params_chen2020NE.c_max - c_n_surf) ** 0.5
# R_ct_n = V_ref / j0_n / (params_chen2020NE.V_ed * params_chen2020NE.epsilon * params_chen2020NE.Sa_vol)

# I = 5
R_ct_n = 2 * V_ref / I * np.arcsinh(I / (params_chen2020NE.V_ed * params_chen2020NE.epsilon * params_chen2020NE.Sa_vol) / 2 / j0_n);

#%%
# calculate terminal voltage
sto = c_p_surf / params_chen2020PE.c_max
OCP_pe = ( 
        -0.8090 * sto + 4.4875 
        - 0.0428 * np.tanh(18.5138 * (sto - 0.5542))
        - 17.7326 * np.tanh(15.7890 * (sto - 0.3117))
        + 17.5842 * np.tanh(15.9308 * (sto - 0.3120))
    )

sto = c_n_surf / params_chen2020NE.c_max
OCP_ne = (
        1.9793 * np.exp(-39.3631 * sto)
        + 0.2482
        - 0.0909 * np.tanh(29.8538 * (sto - 0.1234))
        - 0.04478 * np.tanh(14.9159 * (sto - 0.2769))
        - 0.0205 * np.tanh(30.4444 * (sto - 0.6103))
    )

#%%
# f = open("test\ocp.txt", "w")
# np.savetxt(f, np.c_[sto, OCP_pe, OCP_ne], fmt='%1.6e', delimiter=', ')
# f.close()
I = data['current'].values
V_t = OCP_pe - OCP_ne - I * (R_ct_p + R_ct_n)

#%%
from scipy import integrate
I = data['current'].values
t = data['time'].values
Q = integrate.cumtrapz(I, t, initial=0) / 3600



