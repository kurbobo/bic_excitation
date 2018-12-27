import numpy as np
from numpy.fft import *
from scipy.integrate import odeint
import matplotlib.pyplot as plt
def odeintz(func, z0, t, **kwargs):
    """An odeint-like function for complex valued differential equations."""

    # Disallow Jacobian-related arguments.
    _unsupported_odeint_args = ['Dfun', 'col_deriv', 'ml', 'mu']
    bad_args = [arg for arg in kwargs if arg in _unsupported_odeint_args]
    if len(bad_args) > 0:
        raise ValueError("The odeint argument %r is not supported by "
                         "odeintz." % (bad_args[0],))

    # Make sure z0 is a numpy array of type np.complex128.
    z0 = np.array(z0, dtype=np.complex128, ndmin=1)

    def realfunc(x, t, *args):
        z = x.view(np.complex128)
        dzdt = func(z, t, *args)
        # func might return a python list, so convert its return
        # value to an array with type np.complex128, and then return
        # a np.float64 view of that array.
        return np.asarray(dzdt, dtype=np.complex128).view(np.float64)

    result = odeint(realfunc, z0.view(np.float64), t, **kwargs)

    if kwargs.get('full_output', False):
        z = result[0].view(np.complex128)
        infodict = result[1]
        return z, infodict
    else:
        z = result.view(np.complex128)
        return z
def curcuit(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma, dist, coef,phi,g):
    def model(z, t, Pi):
        U_1 = z[0]
        V_1 = z[1]
        U_2 = z[2]
        V_2 = z[3]
        U_C = z[4]
        V_C = -gamma_C*U_C+sigma*(U_1-U_2)-Pi
        dU_1dt =V_1-gamma_U*U_1+ (chi_2 * U_1 + chi_3 * U_1 ** 2)*V_1-V_C
        dV_1dt =-U_1-gamma_V*V_1
        dU_2dt =V_2-gamma_U*U_2+ (chi_2 * U_2 + chi_3 * U_2 ** 2)*V_2+V_C
        dV_2dt = -U_2-gamma_V*V_2
        dU_Cdt = -gamma_C*U_C+sigma*(U_1-U_2)-Pi
        dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
        return dzdt

    # initial condition
    z_0 = [0.000001, 0, 0.000001, 0, 0]

    # number of time points

    # time points
    # dist = 10000
    # coef = 5
    n = dist * coef
    t = np.linspace(0, dist, n)

    # step input
    Pi = g * np.cos(t+phi) * np.cos(delta * t)  # *np.heaviside(50-t,0)
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
        z = odeint(model, z_0, tspan, args=(Pi[i],))
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




def zfunc(z, t,delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma,g):
    B_1,B_2 = z
    Sigma=sigma*(gamma_C*1j-1)/(1+gamma_C**2)
    gamma=(gamma_V+gamma_U)/2
    alpha=3*chi_3/4+chi_2**2/6
    p=np.sqrt(2)*g*1j
    return [(-1j*delta-gamma + 2j*alpha*( np.abs(B_1)**2 +2*np.abs(B_2)**2))*B_1 +2j*alpha*B_2**2*np.conj(B_1),
            (-1j * delta - gamma +Sigma + 2j * alpha * (np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2))*B_2 +
    2j * alpha * B_1 ** 2 * np.conj(B_2) + p]
# Set up the inputs and call odeintz to solve the system.
z0 = np.array([0.000001, 0])
delta=0.5
chi_2=2e-02*1e02/2
chi_3=1.2e-02*1e02/2
gamma_U=6e-04*1e02
gamma_V=6e-04*1e02
gamma_C=1e-04
sigma=2.5e-04*1e03
dist=250000
coef=1
phi=0
#upper bound - 8.5e-02
t = np.arange(0, dist,1)
gamma=(gamma_V+gamma_U)/2
alpha=3*chi_3/4+chi_2**2/6
g=0.185e-02*1e02/np.sqrt(2)
p=np.sqrt(2)*g#*1j
# Gamma= 1e-03
# alpha= 7.500166666666668e-06
# sigma= 0.0004999995000005001
# sigma2= 4.999995000005e-07
# delta_0=1e-04
print('Gamma=',gamma)
print('alpha=',alpha)
print('sigma=',sigma/2/(1+gamma_C**2))
print("sigma2=", sigma*gamma_C / 2/(1 + gamma_C ** 2))
print('delta_0=',delta)
print('p=', p)
print((sigma*gamma_C /(1 + gamma_C ** 2)+delta)**2-3*(gamma+sigma/(1+gamma_C**2))**2)
# real_curcuit(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma, dist, coef,phi,g):
# t_1,B_1r,B_2r=curcuit(delta=delta,
#                                   chi_2=chi_2,
#                                   chi_3=chi_3,
#                                   gamma_U=gamma_U,
#                            gamma_V=gamma_V,
#                                   gamma_C=gamma_C,
#                            sigma=sigma,
#                                   dist=dist,
#                                   coef=coef,
#                                   phi=phi,
#                            g=g)
z, infodict = odeintz(zfunc, z0, t, args=(delta, chi_2, chi_3, gamma_U,gamma_V,gamma_C,sigma,g), full_output=True)
print(np.abs(z[:,0][2400000]),'=itv')
# if(np.abs(z[:,0][2400000])>=0.25*np.abs(z[:,1][2400000])):
#     print('p=', p)
#     print('g=', g)
#     break
# g-=1e-02
# # import matplotlib.pyplot as plt
# plt.clf()
# max=max(np.abs(B_2r[5000:]))
# lst=[]
# for i in range(dist):
#     lst.append(max)
#print slow amplitudes
plt.plot(t,np.abs(z[:,0]), label='B_1')
plt.plot(t,np.abs(z[:,1]), label='B_2')
# plt.plot(t,lst, label='B_2')
#print real values
# plt.plot(t_1, np.abs(B_1r), label='B_1 real')
# plt.plot(t_1, np.abs(B_2r), label='B_2 real')
# plt.plot(t_1, lst, label='B_2 real')
plt.xlabel('t')
plt.grid(True)
plt.legend(loc='best')
plt.show()
