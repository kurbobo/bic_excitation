"""this scripo is approvement of parameters eqvivalence"""
import numpy as np
import modeling.approx as apx

dist=20000
delta=0.5
chi_2=0
chi_3=0
R=0.10897247358851686
r=0.1
zed=0.1
R_=0.002
k=4#k=C_/C
varepsilon_0=0.142552073818657
coef=1
phi=0


gamma_U=zed/r
gamma_V=R_/zed
gamma_C=k*zed/R
sigma=zed/R
g=zed*varepsilon_0/R

if __name__ == "__main__":
    # Set up the inputs and call odeintz to solve the system.
    z0 = np.array([0.00001, 0])
    t = np.arange(0, dist,1)
    # real_circuit(delta, chi_2, chi_3, R, r,zed, R_, C,  varepsilon_0, dist, coef):
    t_1,B_1r,B_2r=apx.real_circuit(delta=delta,#delta
                                  chi_2=chi_2,
                                   chi_3=chi_3,
                                  R=R,
                                  r=r,
                                  zed=zed,
                                  R_=R_,
                                  k=k,#k=C_/C
                                  varepsilon_0=varepsilon_0,
                                  dist=dist,
                                  coef=coef,
                                  phi=phi)

    t, B_1, B_2 = apx.circuit(delta=delta,
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


    # z, infodict = apx.odeintz(apx.zfunc, z0, t, args=(gamma,sigma,sigma2,alpha,p,delta), full_output=True)
    import matplotlib.pyplot as plt
    plt.clf()
    # #print slow amplitudes
    # plt.plot(t,np.abs(z[:,0]), label='B_1')
    # plt.plot(t,np.abs(z[:,1]), label='B_2')
    # #print real values
    plt.plot(t_1, np.abs(B_1r), label='B_1 real')
    plt.plot(t_1, np.abs(B_2r), label='B_2 real')
    # plt.plot(t, np.abs(B_1), label='B_1 trans')
    # plt.plot(t, np.abs(B_2), label='B_2 trans')
    plt.xlabel('t')
    plt.grid(True)
    plt.legend(loc='best')
    plt.show()
