import numpy as np
import matplotlib.pyplot as plt
import modeling.approx as apx
"""this script is an attempt to prove that we get 
approximately same decisions before and after slow amplitudes approximation"""
# Set up the inputs and call odeintz to solve the system.
z0 = np.array([0.000001, 0])
delta= 0.5
chi_2= 1.0
chi_3= 0.01
sigma= 0.25
gamma_U= 9.999999999999999e-05
gamma_V= 9.999999999999999e-05
gamma_C= 0.05
g= 0.005
dist=250000
phi=0
coef=1
t = np.arange(0, dist,1)
gamma=(gamma_V+gamma_U)/2
alpha=3*chi_3/4+chi_2**2/6
p=np.sqrt(2)*g
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
#print slow amplitudes
plt.clf()
plt.plot(t_1, np.abs(B_1r), label='B_1 real')
plt.plot(t_1, np.abs(B_2r), label='B_2 real')
plt.plot(t,np.abs(z[:,0]), label='B_1')
plt.plot(t,np.abs(z[:,1]), label='B_2')
#print real values
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.savefig('pics/hope.png')
