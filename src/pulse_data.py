
import numpy as np
import pandas as pd
from scipy import interpolate

import codecs

# %%
# --------------------------------------------------------------------------------
# use numpy loadtxt
# --------------------------------------------------------------------------------
# # %% load data 
# import numpy as np
# import codecs

# fname = './icell/1Cpulse_cell1.txt'
# filecp = codecs.open(fname, encoding = 'cp1252')
# data = np.loadtxt(filecp, delimiter='\t', skiprows=91, usecols=[8,10,11])

# #%% plot data
# # fig, ax1 = plt.subplots(tight_layout=True)
# # ax1.plot(data[:,0], data[:,1], 'C0-', label='volt')

# # ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
# # ax2.plot(data[:,0], data[:,2], 'C1')

# # plt.show()


# --------------------------------------------------------------------------------
# use panda 
# https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.read_table.html
# pandas.read_table
# pandas.read_csv
# --------------------------------------------------------------------------------

# import pandas as pd
# file = 'icell/1Cpulse_cell1.txt' 
# df = pd.read_table(file, sep='\t', skiprows=91, header='infer', usecols=[8,10,11] )


class pulseData:
    """
    Data from battery cell pulse charge/discharge test.
    cycle 0: constant current, one discharge pulse 
    cycle 1: each pulse 4% capacity 25 pulses 
    cycle 2: each pulse 20% capacity 5 pulses
    """

    def __init__(self, path):
        """
        Initialize with path to pulse data file.

        Parameters
        ----------
        path : str
            Path to pulse data file. 

        Attributes
        ----------
        time : vector
            Time vector for HPPC battery cell test data [s]
        current : vector
            Current from HPPC battery cell during test [A]
        voltage : vector
            Voltage from HPPC battery cell during test [V]
        cycleNr : vector
            cycleNr indicating for constant current, 4% pulse, 20% pulse [-]
        """

        # loadtxt slow
        # filecp = codecs.open(path, encoding = 'cp1252')
        # data = np.loadtxt(filecp, delimiter='\t', skiprows=91, usecols=[8,10,11,29])
        # # 8---10---11---29
        # # time---volt---curr---cycle number
        # self.time     = data[:,0]
        # self.current  = data[:,2]
        # self.voltage  = data[:,1]
        # self.cycleNrs = data[:,3]


        # # df = pd.read_csv(path,sep='\t',skiprows=91,usecols=[8,10,11,29]) 
        # df = pd.read_table(path,sep='\t',skiprows=91,usecols=[8,10,11,29])
        # # self.time     = df['time/s'].values
        # # self.current  = df['I/mA'].values
        # # self.voltage  = df['Ecell/V'].values
        # # self.cycleNrs = df['cycle number'].values
        # # 8---10---11---29
        # # time---volt---curr---cycle number
        # data = df.values 
        # self.time     = data[:,0]
        # self.current  = data[:,2]
        # self.voltage  = data[:,1]
        # self.cycleNrs = data[:,3]

        # %% 
        tempfram = pd.read_table(path, skiprows=1, nrows=1, header = None)
        secline = tempfram[0][0] 
        
        if secline[0:15] == 'Nb header lines':
            secline_spt = secline.split()
            headlines = int( secline_spt[-1] )
            skipines = headlines - 1
            
            filecp = codecs.open(path, encoding = 'cp1252')
            df = pd.read_table(filecp, sep='\t', skiprows=skipines, header=0)

            self.time     = df['time/s'].values
            self.current  = df['I/mA'].values / 1000 # change to A
            self.voltage  = df['Ecell/V'].values
            self.cycleNrs = df['cycle number'].values
        else:
            print('go back to check header lines!')
        
        #%%

        # number of C rates
        # self.xc = xc



    @classmethod
    def process_one_cycle(cls, path, cycNr):
        """
        retrieve data for one cycle

        cycNr should be 0, 1, 2
        """

        data = cls(path)

        ids = np.where(data.cycleNrs == cycNr)[0]

        # data.time     = data.time[ids] - data.time[ids[0]]
        data.time     = data.time[ids]
        data.current  = data.current[ids]
        data.voltage  = data.voltage[ids]
        data.cycleNrs = data.cycleNrs[ids]

        return data

    @classmethod
    def process_all_cycles(cls, path):
        """
        retrieve data for all cycle

        cycNr should be 0, 1, 2
        """

        data_0 = cls(path) 
        data_1 = cls(path) 
        data_2 = cls(path) 

        ids_0 = np.where(data_0.cycleNrs == 0)[0]
        ids_1 = np.where(data_1.cycleNrs == 1)[0]
        ids_2 = np.where(data_2.cycleNrs == 2)[0]

        data_0.time     = data_0.time[ids_0] - data_0.time[ids_0[0]]
        data_0.current  = data_0.current[ids_0]
        data_0.voltage  = data_0.voltage[ids_0]
        data_0.cycleNrs = data_0.cycleNrs[ids_0]


        data_1.time     = data_1.time[ids_1] - data_1.time[ids_1[0]]
        data_1.current  = data_1.current[ids_1]
        data_1.voltage  = data_1.voltage[ids_1]
        data_1.cycleNrs = data_1.cycleNrs[ids_1]


        data_2.time     = data_2.time[ids_2] - data_2.time[ids_2[0]]
        data_2.current  = data_2.current[ids_2]
        data_2.voltage  = data_2.voltage[ids_2]
        data_2.cycleNrs = data_2.cycleNrs[ids_2]

        return data_0, data_1, data_2

    # @staticmethod
    def get_indices_discharge(self):
        """
        detect the 5 special points of pulse discharge and
        remove the prevailing data
        """
        current  = self.current
        # voltage  = data.voltage

        n = len(current)

        initiPt = []    # first point of each pulse also end point of last point
        discSta = []    # discharge period starting point
        discEnd = []    # discharge period end point
        restSta = []    # rest period starting point
        # restEnd       # rest period end point also next pulse first point

        # minimum current abs
        dI = abs(min(current))

        for k in range(1, n):
            if abs(current[k-1] - current[k] - dI) < dI/100:
                initiPt.append(k-1)
                discSta.append(k)
            if abs(current[k] - current[k-1] - dI) < dI/100:
                discEnd.append(k-1)
                restSta.append(k)

        restEnd = initiPt.copy()

        # locate the last point

        # way 1: simply the last experimental data
        # restEnd.append(n-1) 
        # way 2: make sure the last pulse has same time length as the last second pulse
        # restEnd.append(restSta[-1] + restEnd[-1]-restSta[-2])
        # way 3: make sure the last pulse has same time length as the last second pulse
        # this is for test with 25 pulses
        if len(restEnd) > 1:
            deltaTime = self.time[restEnd[-1]] - self.time[restSta[-2]]
            for i in range(restSta[-1]+1, n):
                if abs(self.time[i] - self.time[restSta[-1]] - deltaTime) < 10.0:
                    restEnd.append(i)
                    break
        # this is for constant current
        else:
            restEnd.append(n-1) 

        # remove the first value
        del restEnd[0]

        initiPt = np.array(initiPt)
        discSta = np.array(discSta)
        discEnd = np.array(discEnd)
        restSta = np.array(restSta)
        restEnd = np.array(restEnd)


        allids = [initiPt, discSta, discEnd, restSta, restEnd]

        return allids

    # remove redundant data before the first pulse
    def trim(self, allids):

        initiPt = allids[0]
        discSta = allids[1] # also could be for charge case
        discEnd = allids[2] # also could be for charge case
        restSta = allids[3]
        restEnd = allids[4]

        self.time      = self.time[initiPt[0]:restEnd[-1]+1] - self.time[initiPt[0]]
        self.current   = self.current[initiPt[0]:restEnd[-1]+1]
        self.voltage   = self.voltage[initiPt[0]:restEnd[-1]+1]
        self.cycleNrs  = self.cycleNrs[initiPt[0]:restEnd[-1]+1] 

        discSta = discSta - initiPt[0]
        discEnd = discEnd - initiPt[0]
        restSta = restSta - initiPt[0]
        restEnd = restEnd - initiPt[0]
        # initiPt goes the last coz used above
        initiPt = initiPt - initiPt[0]

        allids_updated = [initiPt, discSta, discEnd, restSta, restEnd]

        return allids_updated


    def soc(self, params):
        """
        state of charge by Coulomb Counting
        q_cell in Ah
        """
        current = self.current
        time = self.time

        Q_cell = params.Q_cell
        dt = np.diff(time)

        nc = len(current)
        z = np.ones(nc)

        for k in range(1, nc):

            # current step
            # z[k] = z[k-1] + current[k] * dt[k-1] / Q_cell

            # previous step
            # z[k] = z[k-1] - current[k-1] * dt[k-1] / Q_cell

            # previous + current step
            z[k] = z[k - 1] + (current[k] +
                               current[k - 1]) / 2 * dt[k - 1] / Q_cell

        return z


    def get_indices_driveCycle(self):
        """
        detect the 5 special points of pulse discharge and
        remove the prevailing data
        """
        current  = self.current
        # voltage  = data.voltage

        n = len(current)

        initiPt = []    # first point of the drive cycle
        drivSta = []    # discharge period starting point
        drivEnd = []    # discharge period end point
        restSta = []    # rest period starting point
        restEnd = []    # rest period end point 

        # minimum current abs
        dI = abs(max(current))

        for k in range(1, n):
            if current[k-1] - current[k] - dI > 0:
                initiPt.append(k-1)
                drivSta.append(k)
            if current[k] - current[k-1] - dI > 0:
                drivEnd.append(k-1)
                restSta.append(k)

        restEnd.append(n-1) 

        initiPt = np.array(initiPt)
        drivSta = np.array(drivSta)
        drivEnd = np.array(drivEnd)
        restSta = np.array(restSta)
        restEnd = np.array(restEnd)


        allids = [initiPt, drivSta, drivEnd, restSta, restEnd]

        return allids

    def add_1stpoint_charge(self, data_pre):
        """
        add the last point of data_0 as the first point of data_1
        """
        
        self.time      = np.append( data_pre.time[-1],    self.time    )
        self.current   = np.append( data_pre.current[-1], self.current )
        self.voltage   = np.append( data_pre.voltage[-1], self.voltage )
        self.cycleNrs  = np.append( self.cycleNrs[0],     self.cycleNrs)

        return None

    def get_indices_charge(self):
        """
        detect the 5 special points of pulse discharge and
        remove the prevailing data
        """
        current  = self.current
        # voltage  = data.voltage

        n = len(current)

        initiPt = []    # first point of each pulse also end point of last point
        chagSta = []    # discharge period starting point
        chagEnd = []    # discharge period end point
        restSta = []    # rest period starting point
        # restEnd       # rest period end point also next pulse first point

        # minimum current abs
        dI = abs(max(current))

        for k in range(1, n):
            if abs(current[k] - current[k-1] - dI) < dI/100:
                initiPt.append(k-1)
                chagSta.append(k)
            if abs(current[k-1] - current[k] - dI) < dI/100:
                chagEnd.append(k-1)
                restSta.append(k)

        restEnd = initiPt.copy()

        # locate the last point

        # way 1: simply the last experimental data
        # restEnd.append(n-1) 
        # way 2: make sure the last pulse has same time length as the last second pulse
        # restEnd.append(restSta[-1] + restEnd[-1]-restSta[-2])
        # way 3: make sure the last pulse has same time length as the last second pulse
        # this is for test with 20 pulses
        if len(restEnd) > 1:
            deltaTime = self.time[restEnd[-1]] - self.time[restSta[-2]]
            for i in range(restSta[-1]+1, n):
                if abs(self.time[i] - self.time[restSta[-1]] - deltaTime) < 10.0:
                    restEnd.append(i)
                    break
        # this is for constant current
        else:
            jump = abs(min(current))
            for k in range(1, n):
                if abs(current[k-1] - current[k] - jump ) < dI/100:
                    restEnd.append(k-1)

        # remove the first value
        del restEnd[0]

        initiPt = np.array(initiPt)
        chagSta = np.array(chagSta)
        chagEnd = np.array(chagEnd)
        restSta = np.array(restSta)
        restEnd = np.array(restEnd)


        allids = [initiPt, chagSta, chagEnd, restSta, restEnd]

        return allids


    def trim_CCcharge(self, allids):

        initiPt = allids[0]
        chagSta = allids[1]
        chagEnd = allids[2]
        restSta = allids[3]
        restEnd = allids[4]

        self.time      = self.time[initiPt[0]:chagEnd[-1]+1] - self.time[initiPt[0]]
        self.current   = self.current[initiPt[0]:chagEnd[-1]+1]
        self.voltage   = self.voltage[initiPt[0]:chagEnd[-1]+1]
        self.cycleNrs  = self.cycleNrs[initiPt[0]:chagEnd[-1]+1] 

        chagSta = chagSta - initiPt[0]
        chagEnd = chagEnd - initiPt[0]
        restSta = restSta - initiPt[0]
        restEnd = restEnd - initiPt[0]
        # initiPt goes the last coz used above
        initiPt = initiPt - initiPt[0]

        allids_updated = [initiPt, chagSta, chagEnd, restSta, restEnd]

        return allids_updated