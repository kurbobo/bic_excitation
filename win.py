import approx as apx
import numpy as np, matplotlib.pyplot as plt
import os
xi=1#np.sqrt(2)
def zfunc_nonlean_B(z, t, delta, alpha, gamma, gamma_C, s, g):
    b, B_20 = z
    # alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b + 1j * alpha*(np.abs(b) ** 2+2 * np.abs(B_20) ** 2) * b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta- gamma)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + 1j / (1j + gamma_C) * (-s * B_20 + g / 2 / xi)]
z0 = np.array([3.5333586086881496e-16+(1.063600516608418e-16j),0.5508949285073118+(-0.13223625528845565j)])#0.5601364600973774, 0.33891810476589906
delta=0.8
alpha=1
gamma=0.1
gamma_C=1
s=1
#39500
dist=1000
phi=0
gamma_light=gamma + s/(1+gamma_C**2)
delta_light=delta + gamma_C*s/(1+gamma_C**2)
print('delta_light = ',delta_light)
print('gamma_light = ',gamma_light)
g=1.95#15.495031461729917
p=1j*g/2/xi/(1j+gamma_C)
print('p = ',p)
coef=5
start=0*coef#dist-100
stop= dist*coef
t = np.arange(0, dist,1/coef)
# for g in np.arange(19.92034024566258,19.92034024566258-10,-0.1):
z, infodict = apx.odeintz(zfunc_nonlean_B, z0, t, args=(delta, alpha, gamma, gamma_C, s, g), full_output=True)
B_1=z[:,0]*xi#*np.exp(1j*(1+delta)*t)
B_2=z[:,1]*xi#*np.exp(1j*(1+delta)*t)
print("B_1_stat = ",str(np.real(B_1)[(dist-10)*coef])+'+('+str(np.imag(B_1)[(dist-10)*coef])+'j)')
print("B_2_stat = ",str(np.real(B_2)[(dist-10)*coef])+'+('+str(np.imag(B_2)[(dist-10)*coef])+'j)')#0.4362262696148833
# print("B_2_stat_imag = ",str(np.imag(B_2)[(dist-10)*coef]))
plt.clf()
# plt.text(10000, 0.35, r" $\alpha=$"+str(alpha)+r" $ \gamma_C = $"+str(gamma_C)+r" s = "+str(s) + " g =" + str(g)+r" $\delta = $"+str(delta) + r" $B_{10} |{t=0} = $"+str(z0[0]))
plt.plot(t[start:stop],(np.abs((B_2)))[start:stop], label='B_2_stat = '+ str(np.abs(B_2)[(dist-10)*coef]))
plt.plot(t[start:stop], (np.abs((B_1)))[start:stop], label='B_1_stat = ' + str(np.abs(B_1)[(dist-10)*coef]))
plt.text(500,0.2,r"g = " + str(g))
print(max(np.abs(B_1)))
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
# plt.savefig("gisteresis/down/"+str(mx+1)+".png")
plt.show()



