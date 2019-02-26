import modeling.approx as apx
import numpy as np, matplotlib.pyplot as plt
xi=1/np.sqrt(2)
def zfunc_nonlean(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    b, B_20 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b+ 1j * alpha*(np.abs(b)** 2+2* np.abs(B_20) ** 2) *b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta- gamma)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + 1j / (1j + gamma_C) * (-s * B_20 + g / 2 / xi)]
z0 = np.array([0.000001/np.sqrt(2), 0])
delta=3e-02
chi_2=0#3e-01
chi_3=3e-02
gamma_U=0#3e-01
gamma_V=gamma_U
gamma_C=3e-02
s=3e-02
start=0
dist=40000
stop=dist
phi=0
g=0.0185*4
coef=1
t = np.arange(0, dist,1)

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
B_1=z[:,0]*np.exp(1j*delta*t)
B_2=z[:,1]*np.exp(1j*delta*t)
plt.clf()
# plt.plot(t_1[start:stop], np.abs(B_2r[start:stop]), label='B_2 real')
# plt.plot(t[start:stop],np.abs(np.real(2 * B_2*xi))[start:stop], label='B_2_first')
plt.plot(t_1[start:stop], np.abs(B_1r[start:stop]), label='B_1 real')
plt.plot(t[start:stop],np.abs(np.real(2*B_1*xi ))[start:stop], label='B_1_first')
#*np.exp(1j*t)+np.conj(z[:, 1])*np.exp(-1j*t))
# plt.plot(t[start:stop],np.real(2 * z2[:, 0])[start:stop], label='B_1_sec')*np.exp(1j*t)+np.conj(z[:, 0])*np.exp(-1j*t))
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()