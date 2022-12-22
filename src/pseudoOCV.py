# import matplotlib.pyplot as plt

from src import params_chen2020Cell
from src.pulse_data import pulseData

# from scipy import signal
import numpy as np
from scipy import interpolate
# from scipy import array

def get_ocv():
    
    
    file = 'D:/Clouds/OneDrive - Imperial College London/shellECN/Nialls data/Pseudo-OCV/ocv.mpt'
    
    ocv_data = pulseData(file) 

    current  = ocv_data.current
    n = len(current)

    dI = abs(min(current))
    discSta = []

    for k in range(1, n):
        if abs(current[k-1] - current[k] - dI) < dI/100:
            discSta.append(k-1)
            
    # discSta[0]: max -> 0
    # discSta[1]: 0   -> min
    if len(discSta) == 2:
        
        ocv_data.time     = ocv_data.time[discSta[1]:] - ocv_data.time[discSta[1]]
        ocv_data.current  = ocv_data.current[discSta[1]:]
        ocv_data.voltage  = ocv_data.voltage[discSta[1]:]
        ocv_data.cycleNrs = ocv_data.cycleNrs[discSta[1]:]

    #----------------------------------------------------------------
    # way 1
    ocv_soc = ocv_data.soc(params_chen2020Cell)
    ocv_vol = ocv_data.voltage 

    f = interpolate.interp1d(ocv_soc, ocv_vol)
    z_pts = np.linspace(max(ocv_soc), min(ocv_soc), 10001)
    v_pts = f(z_pts)

    #----------------------------------------------------------------
    # way 2
    # z_pts = ocv_data.soc(params)
    # v_pts = ocv_data.voltage 
    #----------------------------------------------------------------

    # extrapolate
    slope = (v_pts[-4] - v_pts[-1]) / (z_pts[-4] - z_pts[-1])

    npts_extra  = np.ceil( min(z_pts) / (z_pts[-2] - z_pts[-1]) ).astype('int')
    # print(npts_extra)
    z_pts_extra = np.linspace(min(z_pts), 0.0, npts_extra+1)
    z_pts_extra = z_pts_extra[1:]

    v_pts_extra = v_pts[-1] - (z_pts[-1] - z_pts_extra) * slope

    z_pts = np.append(z_pts, z_pts_extra)
    v_pts = np.append(v_pts, v_pts_extra)

    return v_pts, z_pts
