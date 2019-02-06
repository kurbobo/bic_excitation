import modeling.approx as apx
import numpy as np, matplotlib.pyplot as plt
def zfunc_pr(z, t, delta, gamma_U, gamma_V, gamma_C, s, g):
    B_1, B_2 = z
    gamma = (gamma_V + gamma_U) / 2
    # alpha = 3 * chi_3 / 4 + chi_2 ** 2 / 6
    return [(-1j * delta - gamma ) * B_1,
            (-1j * delta - gamma)*B_2 +1j/(1j+gamma_C)*(-s* B_2 + g/2)]
def zfunc_nonlean(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g):
    B_1, B_2 = z
    gamma = (gamma_V + gamma_U) / 2
    alpha =0.5*(3/4*chi_3-chi_2**2/3) #0.04*chi_3
    p = g * 1j/(1j+gamma_C) # *np.sqrt(2)
    return [(-1j * delta - gamma + 2j * alpha * (
            np.abs(B_1) ** 2 + 2 * np.abs(B_2) ** 2)) * B_1 + 2j * alpha * B_2 ** 2 * np.conj(B_1),
            (-1j * delta - gamma -1j/(1j+gamma_C)*s + 2j * alpha * (
                    np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2)) * B_2 +
            2j * alpha * B_1 ** 2 * np.conj(B_2) + p]
z0 = np.array([0.000001, 0])
amplitude=1
delta=0.6e-02*amplitude#0.5e-02*amplitude
chi_2=0.1e-02*amplitude
chi_3=0.01e-02*amplitude
gamma_U=1e-05 * amplitude
gamma_V=0.1e-05 * amplitude
gamma_C=0.5e-02 * amplitude
s=2.5e-05*amplitude
dist=250000
phi=0
g=1e-01*amplitude
coef=1
t = np.arange(0, dist,1)
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
plt.plot(t_1, np.abs(B_1r), label='B_1 real')
plt.plot(t_1, np.abs(B_2r)  , label='B_2 real')
plt.plot(t, np.abs(z[:, 0]), label='B_1')
plt.plot(t, np.abs(z[:, 1]), label='B_2')
# print(max(np.abs(B_1r[50000:]))/np.abs(z[:, 1][50000:]))
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()