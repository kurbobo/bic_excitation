import numpy as np
import matplotlib.pyplot as plt
import modeling.approx as apx
"""this script is an attempt to prove that we get 
approximately same decisions before and after slow amplitudes approximation"""
# Set up the inputs and call odeintz to solve the system.
z0 = np.array([0.000001, 0])
amplitude=1e+02
delta=0.6e-02*amplitude#0.5e-02*amplitude
# chi_2=0.1e-01*amplitude
# chi_3=0.01e-02*amplitude
# gamma_U=1e-06 * amplitude
# gamma_V=0.1e-05 * amplitude
# gamma_C=0.5e-03 * amplitude
sigma=2.5e-03*amplitude

# g=0.05e-03*amplitude
# for gamma_U in np.arange(0.1e-05*amplitude, 1e-05*amplitude,0.1e-05 * amplitude):
for gamma_V in np.arange(1e-10 * amplitude, 10e-06 * amplitude,0.1e-05 * amplitude):#1e-10 * amplitude
    for gamma_C in np.arange(0.5e-03 * amplitude, 5e-03 * amplitude,0.1e-04 * amplitude):#0.5e-03 * amplitude
        for gamma_U in np.arange(1e-10 * amplitude, 10e-06 * amplitude,1e-06 * amplitude):#1e-10 * amplitude
            for g in np.arange(0.05e-03 * amplitude, 0.2e-03*amplitude,0.01e-03 * amplitude):
                for chi_2 in np.arange(0.1e-02*amplitude, 1e-01 * amplitude, 0.1e-02 * amplitude):  # 0.1e-02*amplitude
                    for chi_3 in np.arange(0.01e-02 * amplitude, 1e-02 * amplitude,
                                           0.1e-02 * amplitude):  # 0.01e-02*amplitude
                        dist=25000
                        phi=0
                        coef=1
                        t = np.arange(0, dist,1)
                        gamma=(gamma_V+gamma_U)/2
                        alpha=(3*chi_3/4+chi_2**2/6)
                        p=np.sqrt(2)*g
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
                        print('----------------------')
                        print('delta=',delta)
                        print('chi_2=', chi_2)
                        print('chi_3=', chi_3)
                        print('sigma=', sigma)
                        print("gamma_U=", gamma_U)
                        print("gamma_V=", gamma_V)
                        print("gamma_C=", gamma_C)
                        print('g=',g)
                        print('----------------------')
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
                                                   g=2*np.sqrt(2)*g)
                        z, infodict = apx.odeintz(apx.zfunc_previous, z0, t, args=(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma,g), full_output=True)
                        plt.clf()
                        #print real values
                        # plt.plot(t_1, np.abs(B_1r), label='B_1 real')
                        # plt.plot(t_1, np.abs(B_2r)  , label='B_2 real')
                        # print slow amplitudes
                        # plt.plot(t, np.abs(z[:, 0]), label='B_1')
                        # plt.plot(t, np.abs(z[:, 1]), label='B_2')
                        plt.xlabel('t')
                        plt.grid(True)
                        plt.legend(loc='best')
                        plt.savefig('pics/pic'+str(chi_2)+'_'+str(chi_3)+'_'+str(gamma_V)+'_'+str(gamma_U)+'_'+str(gamma_C)+'_'+str(g)+'.png')
