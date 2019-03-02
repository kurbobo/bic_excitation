import modeling.approx as apx
import numpy as np, matplotlib.pyplot as plt
xi=1/np.sqrt(2)
def zfunc_nonlean_A(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    A_1, A_2 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma) * A_1 + 1j * alpha*np.abs(A_1) ** 2 * A_1 -1j/2/ (1j + gamma_C) * (s *(A_1-A_2)+ g / xi),
            (-1j * delta - gamma) * A_2 + 1j * alpha * np.abs(A_2) ** 2 * A_2 + 1j / 2 / (1j + gamma_C) * (
                        s * (A_1 - A_2) + g / xi)]
def zfunc_nonlean_B(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    b, B_20 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha = np.abs(xi)**2*( 3*chi_3/2-chi_2 ** 2/6)
    return [(-1j * delta- gamma ) * b + 1j * alpha*(np.abs(b) ** 2+2 * np.abs(B_20) ** 2) * b
            + 1j * alpha * B_20 ** 2 * np.conj(b),
            (-1j * delta- gamma)* B_20 + 1j * alpha * ((np.abs(B_20) ** 2 + 2 * np.abs(b) ** 2) * B_20
             + b ** 2 * np.conj(B_20)) + 1j / (1j + gamma_C) * (-s * B_20 + g / 2 / xi)]
z0 = np.array([0.000001/np.sqrt(2), 0])
delta=9e-02
chi_2=0#3e-01
chi_3=0#1
gamma_U=3e-01
gamma_V=gamma_U
gamma_C=3e-1
s=0#3e-02
start=400#39500
dist=500
stop=dist
phi=0
g=0.5
coef=10
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
z, infodict = apx.odeintz(zfunc_nonlean_A, z0, t, args=(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g), full_output=True)
C_1=z[:,0]
C_2=z[:,1]
A_1=C_1*np.exp(1j*(1+delta)*t)*xi
A_2=C_2*np.exp(1j*(1+delta)*t)*xi
# B_1=z[:,0]*np.exp(1j*delta*t)
# B_2=z[:,1]*np.exp(1j*delta*t)
plt.clf()
plt.plot(t_1[start:stop], (U_2r[start:stop]), label='A_2 real')
plt.plot(t[start:stop],2*np.real(A_2)[start:stop], label='A_2_first_real_part')
# plt.plot(t_1[start:stop], (U_1r[start:stop]), label='A_1 real')
# plt.plot(t[start:stop],2*np.real(A_1)[start:stop], label='A_1_first_real_part')
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()