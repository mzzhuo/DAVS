import numpy as np
import pandas as pd
from scipy import constants

class shellVoltSource:
    """
    shell voltage source 

    Parameters
    ----------
    data : 
        current as function of time

    Attributes
    ----------
    current : vector
    time : vector
    -------
    """
    def __init__(self, data, params, nr_layer, sign):
        """
        Initialize with HPPC battery cell data and model parameters.
        """
        # load_data: current with time
        # self.time    = load_data.time
        # self.current = load_data.current

        # ----------------------------------------------------------------
        # Given: current/mass flux at particle surface at discrete time
        # at t = 0, current is null I = 0
        # initial state of concentration (homogeneous)
        # ----------------------------------------------------------------

        # # time steps of 1 second
        # j0 = params.j0
        # # lithium flux [mol/m^2/sec]
        # # discharge and rest
        # jk = np.hstack((j0*np.ones(1800), np.zeros(3600)))
        # # discharge and rest, then charge and rest
        # jk = np.hstack((jk, np.negative(jk)))
        # # jk = np.insert(jk, 0, 0)
        # F = constants.physical_constants["Faraday constant"][0]
        # self.current = jk * 4 * np.pi * params.radius ** 2 * F
        # self.time = np.arange(len(jk)) * 1.0
        # self.time_steps = len(self.time)

        # self.current = data['current density'].values * 4 * np.pi * params.radius ** 2
        self.current = data['current'].values / params.N * sign
        self.time = data['time'].values
        self.time_steps = len(self.time)

        # initial concentration
        self.c_init = params.c_init
        self.c_max = params.c_max
        self.k_coeff = params.k_coeff
        self.Ds = params.Ds
        self.radius = params.radius
        self.nr_layer = nr_layer

    def solve_forward_Euler(self):

        # forward method:
        # c_n^k+1 - c_n^k = dt * I_n^k / F

        current = self.current
        time = self.time
        dt = np.diff(time)
        # dt = np.append(dt, dt[-1])
        # time = np.append(time, time[-1] + dt[-1])

        k_coeff = self.k_coeff
        nr_layer = self.nr_layer
        Ds = self.Ds
        radius = self.radius
        c_init = self.c_init
        time_steps = self.time_steps

        F = constants.physical_constants["Faraday constant"][0]

        # initiate concentration
        # [c_1, c_2, c_3, ..., c_n, c_N]
        c_all = np.zeros((time_steps, nr_layer))
        c_all[0, :] = np.ones((1, nr_layer)) * c_init

        # initiate current
        # [I1, I2, ..., In, I_N]
        I_all = np.zeros((time_steps, nr_layer))

        # get the internal resistances, not function of SoC
        # R_diff from 1 to N-1 stored in vector R_diffs
        # R_diffs = np.zeros(nr_layer-1)
        # for i in range(0, nr_layer - 1):
        #     # R_diffs[i] = k_coeff * nr_layer / (4 * np.pi * (i + 1) ** 2 * radius * Ds)
        #     R_diffs[i] = k_coeff * radius / nr_layer / S / Ds
        dR = radius / nr_layer
        S = 4 * np.pi * (np.arange(1, nr_layer) * radius / nr_layer)**2
        R_diffs = k_coeff * dR / S / Ds

        # vol. of ea. shell * radius ** 3
        shell_vols = 4 / 3 * np.pi * ((np.arange(1, nr_layer + 1) * dR)**3 -
                                      (np.arange(0, nr_layer) * dR)**3)

        # loop over time starting from the first time step
        # to calculate variables at the second time step
        # the first timestep hold initial values
        for k in range(0, time_steps - 1):

            # the driving force is concentration gradient at last timestep
            Volt_all = k_coeff * F * c_all[k, :]

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

            # calculate the residuals
            # first SoCs
            # Vol_par = 4 / 3 * np.pi * radius ** 3
            c_all[k + 1, :] = c_all[k, :] - I_all[k, :] * dt[k] / shell_vols / F

        return time, c_all, I_all

    def solve_backward_Euler(self):

        # forward method:
        # c_n^k - c_n^k-1 = dt * I_n^k / F

        current = self.current
        time = self.time
        dt = np.diff(time)

        k_coeff = self.k_coeff
        nr_layer = self.nr_layer
        Ds = self.Ds
        radius = self.radius
        c_max = self.c_max
        c_init = self.c_init
        time_steps = self.time_steps

        F = constants.physical_constants["Faraday constant"][0]

        # initiate concentration
        # [c_1, c_2, c_3, ..., c_n, c_N]
        c_all = np.zeros((time_steps, nr_layer))
        c_all[0, :] = np.ones((1, nr_layer)) * self.c_init

        # initiate current
        # [I1, I2, ..., In, I_N]
        I_all = np.zeros((time_steps, nr_layer))

        # get the internal resistances, not function of SoC
        # R_diff from 1 to N-1 stored in vector R_diffs
        dR = radius / nr_layer
        S = 4 * np.pi * (np.arange(1, nr_layer) * radius / nr_layer)**2
        R_diffs = k_coeff * dR / S / Ds

        # vol. of ea. shell * radius ** 3
        shell_vols = 4 / 3 * np.pi * ((np.arange(1, nr_layer + 1) * dR)**3 -
                                      (np.arange(0, nr_layer) * dR)**3)

        itrNr = 30
        tol = 1e-12
        res = 1.0

        # rat_Res, rat_vol = shell_vol_ratios.get_ratios(nr_layer)

        # loop over time starting from the second time step
        for k in range(1, time_steps):

            # initial guess to be last step values
            # I_all[k,:] = I_all[k-1,:].copy()

            itr = 0
            c_all[k, :] = c_all[k - 1, :].copy()

            while itr < itrNr:  # and corr > tol:
                itr += 1

                Volt_all = k_coeff * F * c_all[k, :]

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

                # calculate the residuals
                res_vec = (c_all[k, :] - c_all[k - 1, :]
                           ) + I_all[k, :] * dt[k - 1] / shell_vols / F
                res_vec = res_vec / c_max
                res = np.linalg.norm(res_vec) / np.sqrt(len(res_vec))

                if res < tol and itr > 1:
                    break

                c_all[k,:] = c_all[k-1,:] - (I_all[k-1,:] + I_all[k,:]) / 2 * dt[k-1] / shell_vols / F
                # c_all[k, :] = c_all[k - 1, :] - I_all[k, :] * dt[k - 1] / shell_vols / F

        return time, c_all, I_all
