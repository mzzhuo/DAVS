# ----------------------------------------------
# run SPM in PyBaMM: 
# https://github.com/pybamm-team/PyBaMM
# ----------------------------------------------

import pybamm
import numpy as np
import matplotlib.pyplot as plt

import pybamm.mz_develop.output_module as outmod

#%%

model = pybamm.lithium_ion.SPM()
param = pybamm.ParameterValues(chemistry=pybamm.parameter_sets.Chen2020)

#%%
experiment = pybamm.Experiment(
    [
     (
      "Discharge at 2.0 C until 2.5 V", 
      "Rest for 1.0 hour",
      "Charge at 2.0 C until 4.2 V", 
      "Rest for 1.0 hour"
      )
    ] * 1,
    period="0.1 minute",
)

#%% GITT
experiment = pybamm.Experiment(
    ["Rest for 0.02 minutes"] + 
    [
     ("Discharge at 0.5 C for 6 minutes", "Rest for 0.8 hour")
    ] * 20,
    period="0.02 minute",
)

#%% CC discharge

experiment = pybamm.Experiment(
    ["Rest for 0.02 minutes"] + 
    [
     ("Discharge at 2.0 C until 2.5 V")
    ],
    period="0.02 minute",
)

#%%


sim = pybamm.Simulation(
    model, experiment=experiment,
    parameter_values=param,
)


#%%
solution = sim.solve(calc_esoh=False)


#%%
output_variables = outmod.output_variables_spm_par
sim.plot(output_variables,time_unit="minutes")  

#%%
time_in_sec = solution["Time [s]"].entries
j_p = solution["X-averaged positive electrode interfacial current density [A.m-2]"].entries

f = open("pybamm\mz_develop\output\j_p_30s.txt", "w")
np.savetxt(f, np.c_[time_in_sec, j_p], fmt='%1.6e', delimiter=', ')
f.close()

#%%
time_in_sec = solution["Time [s]"].entries
c_p_s = solution["X-averaged positive particle surface concentration [mol.m-3]"].entries

f = open("pybamm\mz_develop\output\c_p_s_6s.txt", "w")
np.savetxt(f, np.c_[time_in_sec, c_p_s], fmt='%1.6e', delimiter=', ')
f.close()

#%%
time_in_sec = solution["Time [s]"].entries
eta_p = solution["X-averaged positive electrode reaction overpotential [V]"].entries

f = open("pybamm\mz_develop\output\eta_p_6s.txt", "w")
np.savetxt(f, np.c_[time_in_sec, eta_p], fmt='%1.6e', delimiter=', ')
f.close()


#%%
time_in_sec = solution["Time [s]"].entries
ocp_p = solution["X-averaged positive electrode open circuit potential [V]"].entries

f = open("pybamm\mz_develop\output\OCP_p_6s.txt", "w")
np.savetxt(f, np.c_[time_in_sec, ocp_p], fmt='%1.6e', delimiter=', ')
f.close()

#%%
time_in_sec = solution["Time [s]"].entries
c_p = solution["X-averaged positive particle concentration [mol.m-3]"].entries

r_p = np.linspace(0.025, 0.975, 20)

f = open("pybamm\mz_develop\output\c_p_6s-pb.txt", "w")
# np.savetxt(f, np.c_[time_in_sec, np.transpose(c_p)], fmt='%1.6e', delimiter=', ')
np.savetxt(f, np.c_[r_p, c_p], fmt='%1.6e', delimiter=', ')
f.close()


#%%
time_in_sec = solution["Time [s]"].entries
I = solution["Current [A]"].entries

f = open("pybamm\mz_develop\output\current_CC_3C.txt", "w")
np.savetxt(f, np.c_[time_in_sec, I], fmt='%1.6e', delimiter=', ')
f.close()

#%%
time_in_sec = solution["Time [s]"].entries
OCV = solution["Terminal voltage [V]"].entries

f = open("pybamm\mz_develop\output\OCV_gitt_pb.txt", "w")
np.savetxt(f, np.c_[time_in_sec[0:-1:20], OCV[0:-1:20]], fmt='%1.6e', delimiter=', ')
f.close()


#%%
Vt = solution["Terminal voltage [V]"].entries
Q = solution["Discharge capacity [A.h]"].entries

f = open("pybamm\mz_develop\output\OCV_CC_2C_pb.txt", "w")
np.savetxt(f, np.c_[Q[0:-1:3], Vt[0:-1:3]], fmt='%1.6e', delimiter=', ')
# np.savetxt(f, np.c_[Q, Vt], fmt='%1.6e', delimiter=', ')
f.close()






