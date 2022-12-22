import numpy as np


# particle radius [m]
radius = 5.22E-6
# maximum concentration [mol/m^3]
c_max = 63104
# initial concentration [mol/m^3]
c_init = 17038
# solid diffusivity, [m^2/s]
Ds = 4E-15;
# coefficient [V.m^3/C] = [Ohm.m^3/s]
k_coeff = 1.0e0

# PE volume
V_ed = 75.6E-6 * 6.5E-2 * 1.58 
# PE active material volume fraction
epsilon = 0.665
# PE particle volume 
V_p = 4 / 3 * np.pi * (radius ** 3)

# surface area to volume ratio - particle
Sa_vol = 3 / radius

# Number of particle
# N = 75.6E-6 * 6.5E-2 * 1.58 * 0.665 / (4 / 3 * pi * 5.22E-6**3)
N = 8665901853.69631