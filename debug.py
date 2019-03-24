import approx as apx
import numpy as np, matplotlib.pyplot as plt
xi=1/2#np.sqrt(2)
def zfunc_nonlean_A(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    A_1, A_2 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha = np.abs(xi)**2*( chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma) * A_1 + 1j * alpha*np.abs(A_1) ** 2 * A_1 -1j/2/ (1j + gamma_C) * (s *(A_1-A_2)+ g / xi),
            (-1j * delta - gamma) * A_2 + 1j * alpha * np.abs(A_2) ** 2 * A_2 + 1j / 2 / (1j + gamma_C) * (
                        s * (A_1 - A_2) + g / xi)]
def zfunc_nonlean_B(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    b, B_20 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha = np.abs(xi)**2*(chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b + 1j * alpha*(np.abs(b) ** 2+2 * np.abs(B_20) ** 2) * b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta- gamma)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + 1j / (1j + gamma_C) * (-s * B_20 + g / 2 / xi)]
z0 = np.array([0.000001/np.sqrt(2), 0])
delta=0.1
chi_3 = 0.3
chi_2 = 0.3
gamma_U=6e-04
gamma_V=gamma_U
gamma_C=1
s=6e-03
#39500
coef=10
dist=5000*coef
start=(dist-100)*coef#dist-100
stop= dist*coef#dist
phi=0
g=0.0185
t = np.arange(0, dist,1/coef)
t_1,U_1r,U_2r=apx.circuit(delta=delta,
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
z, infodict = apx.odeintz(zfunc_nonlean_B, z0, t, args=(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g), full_output=True)
# z2, infodict2 = apx.odeintz(zfunc_nonlean_B2, z0, t, args=(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g), full_output=True)
B_1=z[:,0]*xi#*np.exp(1j*(1+delta)*t)
B_2=z[:,1]*xi#*np.exp(1j*(1+delta)*t)
plt.clf()
plt.plot(t_1[start:stop], (U_2r[start:stop]), label='A_2 real')
# plt.plot(t[start:stop], (2 * np.real(B_1 + B_2))[start:stop], label='B_1_plus_B_2_first')
plt.plot(t_1[start:stop], (U_1r[start:stop]), label='A_1 real')
# plt.plot(t[start:stop],(2*np.real(B_1-B_2))[start:stop], label='B_1_minus_B_2_first')
plt.plot(t[start:stop], (2 * np.abs(B_1 + B_2))[start:stop], label='B_1_plus_B_2_first')
plt.plot(t[start:stop],(2*np.abs(B_1-B_2))[start:stop], label='B_1_minus_B_2_first')
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.annotate(s="delta="+str(delta)+"chi3="+str(chi_3)+"gamma="+str(gamma_U)+"gamma_C="+str(gamma_C)+"s="+str(s)+"g="+str(g),xy=(100,0.1))
plt.show()
# plt.savefig("approvement_pics/delta="+str(delta)+"chi3="+str(chi_3)+"gamma="+str(gamma_U)+"gamma_C="+str(gamma_C)+"s="+str(s)+"g="+str(g)+".png")