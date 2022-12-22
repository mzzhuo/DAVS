import numpy as np

def get_ratios(nr_layer):
    """
    nr_layer: number of OCVs, shell/layer number
    nr_res: number of resistance
    nr_res = nr_layer - 1
    mode: how to distribute capacity

    see Zhuo's paper

    """

    # assume the total radius is 1.0
    r_tot = 1.0

    nr_res = nr_layer - 1

    # ----------------------------------------------------
    # mode: equal shell thickness
    # ----------------------------------------------------

    # calculate the outer radius of layer i=1,2,3,...,N
    # r[i-1] = the ith layer external radius 
    r = np.zeros(nr_layer)
    for i in range(0, nr_layer):
        r[i] = r_tot / nr_layer * (i + 1)

    # ref is the 1st layer resistance--R_d,1
    # R_d,2, ..., R_d,nr_layer are defined as their ratios to R_d,1
    rat_res = np.ones(nr_res) 
    for i in range(0, nr_res): 
        rat_res[i] = 1.0 / (i + 1) ** 2

    rat_vol = np.ones(nr_layer) 
    for i in range(1, nr_layer):
        rat_vol[i] = (r[i] ** 3 - r[i-1] ** 3) / (r_tot ** 3)
    rat_vol[0] = (r[0] ** 3) / (r_tot ** 3)

    return rat_res, rat_vol
