import modeling.approx as apx
import numpy as np, matplotlib.pyplot as plt
def zfunc_pr(z, t, delta, gamma_U, gamma_V, gamma_C, s, g):
    B_1, B_2 = z
    gamma = (gamma_V + gamma_U) / 2
    return [(-1j * delta - gamma)* B_1,
            (-1j * delta - gamma)*B_2 + 1j/(1j+gamma_C)*(- s * B_2 + g/np.sqrt(2)*np.heaviside(5000-t,0))]
def zfunc_nonlean(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    B_1 , B_2 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha = 3*chi_3/np.sqrt(2)-chi_2 ** 2/12#3*chi_3/16#/2*np.sqrt(2)+2*chi_2 ** 2/3
    p =g*1j/(1j+gamma_C)/np.sqrt(2)#*np.heaviside(5000-t,0)#g * 1j/(1j+gamma_C)/np.sqrt(2)
    # print(1j/(1j+gamma_C)*s)
    return [(-1j * delta- gamma + 1j * alpha * (np.abs(B_1) ** 2 + 2 * np.abs(B_2) ** 2)) * B_1
            + 1j * alpha * B_2 ** 2 * np.conj(B_1),
            (-1j * delta - gamma + 1j * alpha * (np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2)) * B_2
            + 1j * alpha * B_1 ** 2 * np.conj(B_2) + p - 1j / (1j+gamma_C) * s * B_2]
z0 = np.array([0.000001/np.sqrt(2), 0])
delta=0.002#
chi_2=0#0.6
chi_3=0.4
gamma_U=0#6e-03
gamma_V=gamma_U
gamma_C=0#1e-04
s=0.24e-02
start=0#15000
dist=10000
stop=100
phi=0
g=0.0185*2
coef=1
t = np.arange(0, dist,1)
# for delta in np.arange(0.01,0.1,0.01):
#     for c hi_3 in np.arange(0.01,0.1,0.01):
#         for gamma_V in np.arange(0.01,0.62,0.1):
#             for gamma_C in np.arange(0.01,0.6,0.1):
#                 for s in np.arange(0.01,0.62,0.1):
#                     for g in np.arange(0.01,0.62,0.1):
#                         gamma_U=gamma_V



# p=4*np.sqrt(2)*g
t_1,B_1r,B_2r=apx.circuit(delta=delta,
                                                          chi_2=chi_2,
                                                          chi_3=chi_3,
                                                          gamma_U=gamma_U,
                                                   gamma_V=gamma_V,
                                                          gamma_C=gamma_C,
                                                   sigma=s,
                                                          dist=dist,
                                                          coef=coef,
                                                          phi=phi,
                                                   g=g)
z, infodict = apx.odeintz(zfunc_nonlean, z0, t, args=(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g), full_output=True)
# z, infodict = apx.odeintz(zfunc_pr, z0, t, args=(delta, gamma_U, gamma_V, gamma_C, s, g), full_output=True)
plt.clf()
plt.plot(t_1[start:stop], np.abs(B_2r)[start:stop], label='B_2 real')
plt.plot(t_1[start:stop], np.abs(B_1r)[start:stop], label='B_1 real')
plt.plot(t[start:stop],np.sqrt(2)* np.abs(z[:, 0])[start:stop], label='B_1')
plt.plot(t[start:stop],np.sqrt(2)*np.abs(z[:, 1])[start:stop], label='B_2')
# print(max(np.abs(B_1r[50000:]))/np.abs(z[:, 1][50000:]))
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()
# plt.savefig('pics/pic'+'_'+str(gamma_V)+'_'+str(delta)+'_'+str(s)+'_'+str(chi_3)+'_'+str(gamma_C)+'_'+str(g)+'.png')
# if(np.abs(z[:, 0])[2500]>0.1*np.abs(z[:, 1])[2500]):
#     print(str(gamma_V)+'_'+str(delta)+'_'+str(s)+'_'+str(chi_3)+'_'+str(gamma_C)+'_'+str(g))