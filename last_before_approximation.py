import numpy as np
import matplotlib.pyplot as plt
import modeling.approx as apx
"""this script is an attempt to prove that we get 
approximately same decisions before and after slow amplitudes approximation"""
# Set up the inputs and call odeintz to solve the system.
z0 = np.array([0.000001, 0])
amplitude=1e+02
delta=0.5e-02*amplitude
chi_2=0.5e-02*amplitude#*1e02
chi_3=0.62e-02*amplitude#*1e02
gamma_U=5e-04*amplitude#*1e02
gamma_V=5e-04*amplitude#*1e02
gamma_C=1e-04#????????
sigma=2.5e-03*amplitude#*1e02
dist=250000
phi=0
coef=1
#upper bound - 8.5e-02
t = np.arange(0, dist,1)
gamma=(gamma_V+gamma_U)/2
alpha=3*chi_3/4+chi_2**2/6
g=0.185e-02/np.sqrt(2)*amplitude
p=np.sqrt(2)*g#*1j
# delta=0.5;
# sigma=0.12
#  gamma=0.06;
# alpha=0.5
# sigma2=0
# p=0.185
print('Gamma=',gamma)
print('alpha=',alpha)
print('sigma=',sigma/2/(1+gamma_C**2))
print("sigma2=", sigma*gamma_C / 2/(1 + gamma_C ** 2))
print('delta_0=',delta)
print('p=', p)
print((sigma*gamma_C /(1 + gamma_C ** 2)+delta)**2-3*(gamma+sigma/(1+gamma_C**2))**2)
# real_circuit(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma, dist, coef,phi,g):
t_1,B_1r,B_2r=apx.circuit(delta=delta,
                                  chi_2=chi_2,
                                  chi_3=chi_3,
                                  gamma_U=gamma_U,
                           gamma_V=gamma_V,
                                  gamma_C=gamma_C,
                           sigma=sigma,
                                  dist=dist,
                                  coef=coef,
                                  phi=phi,
                           g=g)
z, infodict = apx.odeintz(apx.zfunc_previous, z0, t, args=(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma,g), full_output=True)
# print(np.abs(z[:,0][240000]),'=itv')
# if(np.abs(z[:,0][2400000])>=0.25*np.abs(z[:,1][2400000])):
#     print('p=', p)
#     print('g=', g)
#     break
# g-=1e-02
# # import matplotlib.pyplot as plt
# plt.clf()
# max=max(np.abs(B_2r[5000:]))
# lst=[]
# for i in range(dist):
#     lst.append(max)
#print slow amplitudes
# plt.plot(t,np.abs(z[:,0]), label='B_1')
# plt.plot(t,np.abs(z[:,1]), label='B_2')
# plt.plot(t,lst, label='B_2')
#print real values
plt.plot(t_1, np.abs(B_1r), label='B_1 real')
plt.plot(t_1, np.abs(B_2r), label='B_2 real')
# plt.plot(t_1, lst, label='B_2 real')
# plt.xlabel('t')
# plt.grid(True)
# plt.legend(loc='best')
# plt.show()
