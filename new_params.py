import approx as apx
import numpy as np, matplotlib.pyplot as plt
import os
xi=1#np.sqrt(2)
def zfunc_nonlean_B(z, t, delta, alpha, gamma, p,delta_light,gamma_light):
    b, B_20 = z
    # alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b + 1j * alpha*(np.abs(b) ** 2+2 * np.abs(B_20) ** 2) * b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta_light- gamma_light)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + p]#*np.heaviside(814-t,0)
z0 = np.array([6.527729810487234e-15+2.197688527407939e-15j,0.5571070388267081-0.13189796679776122j])#0.5601364600973774, 0.33891810476589906
delta=0.8#0.2581988897471611
alpha=1
delta_light=1.3
gamma_light=0.6
gamma=0.1
dist=1000
p=0.4625+0.4625j
coef=5
start=0*coef#dist-100
stop= dist*coef
t = np.arange(0, dist,1/coef)
z, infodict = apx.odeintz(zfunc_nonlean_B, z0, t, args=(delta, alpha, gamma,p, delta_light, gamma_light,), full_output=True)
B_1=z[:,0]*xi
B_2=z[:,1]*xi
print("B_1_stat = ",str(np.abs(B_1)[(dist-10)*coef]))
print("B_2_stat = ",str(np.abs(B_2)[(dist-10)*coef]))
plt.clf()
plt.plot(t[start:stop],(np.abs((B_2)))[start:stop], label='B_2_stat = '+ str(np.abs(B_2)[(dist-10)*coef]))
plt.plot(t[start:stop], (np.abs((B_1)))[start:stop], label='B_1_stat = ' + str(np.abs(B_1)[(dist-10)*coef]))
plt.text(500,0.2,r"p = " + str(p))
print(max(np.abs(B_1)))
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()



