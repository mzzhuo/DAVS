import numpy as np


# particle radius [m]
radius = 5.86E-6
# maximum concentration [mol/m^3]
c_max = 33133
# initial concentration [mol/m^3]
c_init = 29866
# solid diffusivity, [m^2/s]
Ds = 3.3E-14;
# coefficient [V.m^3/C] = [Ohm.m^3/s]
k_coeff = 1.0e0

# PE volume
V_ed = 85.2E-6 * 6.5E-2 * 1.58 
# PE active material volume fraction
epsilon = 0.75
# PE particle volume 
V_p = 4 / 3 * np.pi * (radius ** 3)

# surface area to volume ratio - particle
Sa_vol = 3 / radius

# Number of particle
# N = V_ne * epsilon / V_p
N = 7785559792.48589