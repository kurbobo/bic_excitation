import numpy as np
import matplotlib.pyplot as plt
import modeling.approx as apx,evg_debug as evg
"""this script is an attempt to prove that we get 
approximately same decisions before and after slow amplitudes approximation"""
# Set up the inputs and call odeintz to solve the system.
z0 = np.array([0.000001, 0])

delta=0.5e-02#0.5e-02
# chi_2=0.1e-01*amplitude
# chi_3=0.01e-02*amplitude
gamma_U=1e-06
# gamma_V=0.1e-05 * amplitude
# gamma_C=0.5e-03 * amplitude
s=2.5e-03
chi_3=0
chi_2=0
g=0.05e-03

# for gamma_U in np.arange(0.1e-05*amplitude, 1e-05*amplitude,0.1e-05 * amplitude):
for gamma_V in np.arange(1e-06, 10e-06,0.1e-05):#1e-10 * amplitude
    for gamma_C in np.arange(0.5e-03, 5e-03,0.1e-03):#0.5e-03 * amplitude
        # for s in np.arange(2.5e-03 * amplitude, 25e-03 * amplitude,2.5e-03 * amplitude):#1e-10 * amplitude
            # for g in np.arange(0.5e-02*amplitude, 2e-02*amplitude,0.2e-02 * amplitude):
                for chi_2 in np.arange(0, 1e-01, 0.1e-02):  # 0.1e-02*amplitude
                    for chi_3 in np.arange(0.01e-02, 1e-02,
                                           0.1e-02):  # 0.01e-02*amplitude
                                    dist=250000
                                    phi=0
                                    coef=1
                                    t = np.arange(0, dist,1)
                                    gamma=(gamma_V+gamma_U)/2
                                    # alpha=(3*chi_3/4+chi_2**2/6)
                                    p=g#np.sqrt(2)*
                                    # delta=0.5;
                                    # sigma=0.12
                                    #  gamma=0.06;
                                    # alpha=0.5
                                    # sigma2=0
                                    # p=0.185
                                    print('Gamma=',gamma)
                                    # print('alpha=',alpha)
                                    print('sigma=',s/2/(1+gamma_C**2))
                                    print("sigma2=", s*gamma_C / 2/(1 + gamma_C ** 2))
                                    print('delta_0=',delta)
                                    print('p=', p)
                                    print('----------------------')
                                    print('delta=',delta)
                                    print('chi_2=', chi_2)
                                    # print('chi_3=', chi_3)
                                    print('sigma=', s)
                                    print("gamma_U=", gamma_U)
                                    print("gamma_V=", gamma_V)
                                    print("gamma_C=", gamma_C)
                                    print('g=',g)
                                    print('----------------------')
                                    print((s*gamma_C /(1 + gamma_C ** 2)+delta)**2-3*(gamma+s/(1+gamma_C**2))**2)
                                    # real_circuit(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma, dist, coef,phi,g):

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
                                    z, infodict = apx.odeintz(evg.zfunc_nonlean, z0, t, args=(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, s, g), full_output=True)
                                    plt.clf()
                                    #print real values
                                    plt.plot(t_1, np.abs(B_1r), label='B_1 real')
                                    plt.plot(t_1, np.abs(B_2r)  , label='B_2 real')
                                    # print slow amplitudes
                                    plt.plot(t, np.abs(z[:, 0]), label='B_1')
                                    plt.plot(t, np.abs(z[:, 1]), label='B_2')
                                    print('max for B_1 real is',max(np.abs(B_2r)))
                                    print('max for B_1 is',max(np.abs(z[:, 1])))
                                    print('ratio is',max(np.abs(z[:, 1]))/max(np.abs(B_2r)))
                                    plt.xlabel('t')
                                    plt.grid(True)
                                    plt.legend(loc='best')
                                    plt.savefig('pics/pic'+'_'+str(gamma_V)+'_'+str(gamma_U)+'_'+str(gamma_C)+'_'+str(g)+'.png')
