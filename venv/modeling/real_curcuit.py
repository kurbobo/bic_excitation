import numpy as np
from numpy.fft import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
#
# delta = 0.5
# chi_21 = 0  # 0.21
# chi_31 = 0  # 0.004
# chi_22 = chi_21
# chi_32 = chi_31
# R = 10
# r_1 = 500
# r_2 = r_1
# R_1 = 0.001
# R_2 = R_1
# C_1 = 2E-5
# C_2 = C_1
# C = 1.1E-5
# L_1 = 0.5E-4
# L_2 = L_1
# varepsilon_0 = 0.3140916808149406


def real_curcuit(delta, chi_2, chi_3, R, r,zed, R_, k,  varepsilon_0, dist, coef,phi):
    def model(z, t, varepsilon):
        # zed = np.sqrt(L / C_)
        U_1 = z[0]
        V_1 = z[1]
        U_2 = z[2]
        V_2 = z[3]
        U_C = z[4]
        V_C = zed / R * (U_1 - U_2 - U_C - varepsilon)
        dU_1dt = (1 + chi_2 * U_1 + chi_3 * U_1 ** 2) * (-zed / r * U_1 + V_1 - V_C)
        dV_1dt = (-R_ / zed * V_1 - U_1)
        dU_2dt = ((1 + chi_2 * U_2 + chi_3 * U_2 ** 2) * (-zed / r * U_2 + V_C - V_2))
        dV_2dt = (-R_ / zed * V_2 + U_2)
        dU_Cdt = zed * k / R * (-U_C + (U_1 - U_2) - varepsilon)
        dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
        return dzdt

    # initial condition
    z_0 = [0.00001, 0, 0.00001, 0, 0]

    # number of time points

    # time points
    # dist = 10000
    # coef = 5
    n = dist * coef
    t = np.linspace(0, dist, n)

    # step input
    varepsilon = varepsilon_0 * np.cos(t+phi) * np.cos(delta * t)  # *np.heaviside(50-t,0)
    # store solution
    U_1 = np.empty_like(t)
    V_1 = np.empty_like(t)
    U_2 = np.empty_like(t)
    V_2 = np.empty_like(t)
    U_C = np.empty_like(t)
    # record initial conditions
    U_1[0] = z_0[0]
    V_1[0] = z_0[1]
    U_2[0] = z_0[2]
    V_2[0] = z_0[3]
    U_C[0] = z_0[4]
    # solve ODE
    for i in range(1, n):
        # span for next time step
        tspan = [t[i - 1], t[i]]
        # solve for next step
        z = odeint(model, z_0, tspan, args=(varepsilon[i],))
        # store solution for plotting
        U_1[i] = z[1][0]
        V_1[i] = z[1][1]
        U_2[i] = z[1][2]
        V_2[i] = z[1][3]
        U_C[i] = z[1][4]
        # next initial condition
        z_0 = z[1]
    B_1 = (U_1 + U_2) / 2
    B_2 = (U_1 - U_2) / 2
    return t, B_1, B_2

# # plot results
# plt.plot(t[0:500],np.abs((B_1[:500])),'g',label='dark mode')
# # plt.plot(t[:],varepsilon[:],'r:',label='vareps',)
# plt.plot(t[:500],np.abs((B_2[:500])),'b-',label='light mode')
# # plt.plot(t[:],B_one,'g',label='dark mode')
# # plt.plot(t[:],B_two,'b-',label='light mode')
# plt.ylabel('values')
# plt.text(4000,0.1,'omega= '+str(c))
# plt.xlabel('time')
# plt.legend(loc='best')
# plt.show()
#
# # Gamma= 0.060000000000000005
# alpha= 0.5
# p= 0.18508029784877095
# sigma= 0.13236581096849476
# sigma2= 0.10027712952158695
