#%%
f = open("test/OCV_gitt_ecm.txt", "w")
np.savetxt(f, np.c_[time[0:-1:20], V_t[0:-1:20]], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("test/OCV_CC_3C_ecm.txt", "w")
np.savetxt(f, np.c_[Q[0:-1:4], V_t[0:-1:4]], fmt='%1.6e', delimiter=', ')
# np.savetxt(f, np.c_[Q, V_t], fmt='%1.6e', delimiter=', ')
f.close()

#%%
f = open("test/Rct_pe.txt", "w")
np.savetxt(f, np.c_[sto, c_p_surf, R_ct_p_appro, R_ct_p_5A, R_ct_p_10A], fmt='%1.6e', delimiter=', ')
f.close()


#%%
f = open("test/Rct_ne.txt", "w")
np.savetxt(f, np.c_[sto, c_n_surf, R_ct_n_appro, R_ct_n_5A, R_ct_n_10A], fmt='%1.6e', delimiter=', ')
f.close()