import matplotlib.pyplot as plt

#%%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(time, c_all_ne[:,0]/1000, 'b-')
ax1.plot(time, c_all_ne[:,10]/1000, 'r-')
ax1.plot(time, c_all_ne[:,-1]/1000, 'r-')
ax1.legend()
plt.show()



#%%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(time, j0_p, 'b-')
# ax1.plot(time, j0_n, 'r-')
ax1.legend()
plt.show()



#%%
fig, ax1 = plt.subplots(figsize=(8, 5),tight_layout=True)
ax1.plot(time, V_t, 'b-')
ax1.legend()
plt.show()

