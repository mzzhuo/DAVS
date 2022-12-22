# ----------------------------------------------
# calculate charge transfer resistance
# produce results in Section 3.2
# ----------------------------------------------

import numpy as np
from src import params_chen2020PE
from src import params_chen2020NE
from scipy import constants

R_cons = constants.R
F_cons = constants.physical_constants["Faraday constant"][0]
T_cons = 298.15
V_ref = R_cons * T_cons / F_cons


# calculate particle surface concentration
sto = np.linspace(17038/63104, 0.97, 100)
c_p_surf = 63104 * sto

# calculate exchange current density
j0_p = 3.42e-6 * 1000 ** 0.5 * c_p_surf ** 0.5 * (params_chen2020PE.c_max - c_p_surf) ** 0.5
# calculate the charge transfer resistance
# linear approximation
#%%
R_ct_p_appro = V_ref / j0_p / (params_chen2020PE.V_ed * params_chen2020PE.epsilon * params_chen2020PE.Sa_vol)  

# related to applied current
I = 5
R_ct_p_5A = 2 * V_ref / I * np.arcsinh(I / (params_chen2020PE.V_ed * params_chen2020PE.epsilon * params_chen2020PE.Sa_vol) / 2 / j0_p);

I = 10
R_ct_p_10A = 2 * V_ref / I * np.arcsinh(I / (params_chen2020PE.V_ed * params_chen2020PE.epsilon * params_chen2020PE.Sa_vol) / 2 / j0_p);

#%%
sto = np.linspace(0.04, 29866/33133, 100)
c_n_surf = 33133 * sto
j0_n = 6.48e-7 * 1000 ** 0.5 * c_n_surf ** 0.5 * (params_chen2020NE.c_max - c_n_surf) ** 0.5

R_ct_n_appro = V_ref / j0_n / (params_chen2020NE.V_ed * params_chen2020NE.epsilon * params_chen2020NE.Sa_vol)

I = 5
R_ct_n_5A = 2 * V_ref / I * np.arcsinh(I / (params_chen2020NE.V_ed * params_chen2020NE.epsilon * params_chen2020NE.Sa_vol) / 2 / j0_n);

I = 10
R_ct_n_10A = 2 * V_ref / I * np.arcsinh(I / (params_chen2020NE.V_ed * params_chen2020NE.epsilon * params_chen2020NE.Sa_vol) / 2 / j0_n);