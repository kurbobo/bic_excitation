import modeling.approx as apx
import numpy as np, matplotlib.pyplot as plt
xi=1#np.sqrt(2)
def zfunc_nonlean_B(z, t, delta, alpha, gamma, gamma_C, s, g):
    b, B_20 = z
    # alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b + 1j * alpha*(np.abs(b) ** 2+2 * np.abs(B_20) ** 2) * b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta- gamma)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + 1j / (1j + gamma_C) * (-s * B_20 + g / 2 / xi)]
z0 = np.array([0.00001, 0])
delta=0.2#0.6
alpha=1
gamma=6e-02
gamma_C=1e-01
s=2.59736e+01
#39500
dist=100000
phi=0
g=14.4252
coef=5
start=0*coef#dist-100
stop= dist*coef
t = np.arange(0, dist,1/coef)
# for g in np.arange(19.92034024566258,19.92034024566258-10,-0.1):
z, infodict = apx.odeintz(zfunc_nonlean_B, z0, t, args=(delta, alpha, gamma, gamma_C, s, g), full_output=True)
B_1=z[:,0]*xi#*np.exp(1j*(1+delta)*t)
B_2=z[:,1]*xi#*np.exp(1j*(1+delta)*t)
print("B_1_stat = ",str(np.abs(B_1)[(dist-10)*coef]))
print("B_2_stat = ",str(np.abs(B_2)[(dist-10)*coef]))#0.4362262696148833
plt.clf()
# plt.text(10000, 0.35, r" $\alpha=$"+str(alpha)+r" $ \gamma_C = $"+str(gamma_C)+r" s = "+str(s) + " g =" + str(g)+r" $\delta = $"+str(delta) + r" $B_{10} |{t=0} = $"+str(z0[0]))
plt.plot(t[start:stop],(np.abs((B_2)))[start:stop], label='B_2r')
plt.plot(t[start:stop], (np.abs((B_1)))[start:stop], label='B_1r')
# plt.plot(t[start:stop],np.abs(B_2)[start:stop], label='B_2')
# plt.plot(t[start:stop], np.abs(B_1)[start:stop], label='B_1')
print(max(np.abs(B_1)))
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()
# plt.savefig("pics/delta="+str(delta)+"alpha="+str(alpha)+"gamma="+str(gamma)+"gamma_C="+str(gamma_C)+"s="+str(s)+"g="+str(g)+" B_1_0="+str(z0[0])+" B_2_0="+str(z0[1])+".png")