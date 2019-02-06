import numpy as np
from scipy.integrate import odeint


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


# Define the right-hand-side of the differential equation.
# here dark mode is B_1
def zfunc(z, t, gamma, sigma, sigma2, alpha, p, delta):
    B_1, B_2 = z
    return [(-1j * delta - gamma + 2j * alpha * (
            np.abs(B_1) ** 2 + 2 * np.abs(B_2) ** 2)) * B_1 + 2j * alpha * B_2 ** 2 * np.conj(B_1),
            (-1j * delta - gamma - 2 * sigma + 2j * sigma2 + 2j * alpha * (
                    np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2)) * B_2 +
            2j * alpha * B_1 ** 2 * np.conj(B_2) + p]


# def real_circuit(delta, chi_2, chi_3, R, r, zed, R_, k, varepsilon_0, dist, coef, phi):
#     def model(z, t, varepsilon):
#         # zed = np.sqrt(L / C_)
#         U_1 = z[0]
#         V_1 = z[1]
#         U_2 = z[2]
#         V_2 = z[3]
#         U_C = z[4]
#         V_C = zed / R * (U_1 - U_2 - U_C - varepsilon)
#         dU_1dt = (1 + chi_2 * U_1 + chi_3 * U_1 ** 2) * (-zed / r * U_1 + V_1 - V_C)
#         dV_1dt = (-R_ / zed * V_1 - U_1)
#         dU_2dt = ((1 + chi_2 * U_2 + chi_3 * U_2 ** 2) * (-zed / r * U_2 + V_C - V_2))
#         dV_2dt = (-R_ / zed * V_2 + U_2)
#         dU_Cdt = zed * k / R * (-U_C + (U_1 - U_2) + varepsilon)
#         dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
#         return dzdt
#
#     # initial condition
#     z_0 = [0.00001, 0, 0.00001, 0, 0]
#
#     # number of time points
#
#     # time points
#     # dist = 10000
#     # coef = 5
#     n = dist * coef
#     t = np.linspace(0, dist, n)
#
#     # step input
#     varepsilon = varepsilon_0 * np.cos(t + phi) * np.cos(delta * t)/np.sqrt(2)  # *np.heaviside(100000-t,0)
#     # store solution
#     U_1 = np.empty_like(t)
#     V_1 = np.empty_like(t)
#     U_2 = np.empty_like(t)
#     V_2 = np.empty_like(t)
#     U_C = np.empty_like(t)
#     # record initial conditions
#     U_1[0] = z_0[0]
#     V_1[0] = z_0[1]
#     U_2[0] = z_0[2]
#     V_2[0] = z_0[3]
#     U_C[0] = z_0[4]
#     # solve ODE
#     for i in range(1, n):
#         # span for next time step
#         tspan = [t[i - 1], t[i]]
#         # solve for next step
#         z = odeint(model, z_0, tspan, args=(varepsilon[i],))
#         # store solution for plotting
#         U_1[i] = z[1][0]
#         V_1[i] = z[1][1]
#         U_2[i] = z[1][2]
#         V_2[i] = z[1][3]
#         U_C[i] = z[1][4]
#         # next initial condition
#         z_0 = z[1]
#     B_1 = (U_1 + U_2) / 2
#     B_2 = (U_2 - U_1) / 2
#     return t, B_1, B_2


""""circuit is the same as realcircuit, but earlier variant:
    parameters from realcircuit are here functions of 
    gamma_U, gamma_V, gamma_C and so on"""


def circuit(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, sigma, dist, coef, phi, g):
    def model(z, t, Pi):
        U_1 = z[0]
        V_1 = z[1]
        U_2 = z[2]
        V_2 = z[3]
        U_C = z[4]
        V_C = -gamma_C * U_C + sigma * (U_1 - U_2) + Pi
        dU_1dt = V_1 - gamma_U * U_1 + (chi_2 * U_1 + chi_3 * U_1 ** 2) * V_1 - V_C
        dV_1dt = -U_1 - gamma_V * V_1
        dU_2dt = V_2 - gamma_U * U_2 + (chi_2 * U_2 + chi_3 * U_2 ** 2) * V_2 + V_C
        dV_2dt = -U_2 - gamma_V * V_2
        dU_Cdt = -gamma_C * U_C + sigma * (U_1 - U_2) + Pi
        dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
        return dzdt

    # initial condition
    z_0 = [0.000001, 0, 0.000001, 0, 0]

    # number of time points

    # time points
    # dist = 10000
    # coef = 5
    n = dist * coef  # было умножить
    t = np.linspace(0, dist, n)

    # step input
    Pi =2*g * np.cos((1 + delta) * t + phi)  *np.heaviside(5000-t,0)
    # (np.exp((1j+delta)*t)).real
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


""""zfunc_previous is from last_before_approximation,
    the same as zfunc but there are not sigmas and gammas explicitly,
    but they are as functions of gamma_U,gamma_V and so on"""


def zfunc_previous(z, t, delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, sigma, g):
    B_1, B_2 = z
    i_Sigma2_minus_sigma = sigma * (gamma_C * 1j - 1) /2/ (1 + gamma_C ** 2)
    gamma = (gamma_V + gamma_U) / 2
    alpha = chi_3 /25 #+ chi_2 ** 2 / 6
    p = 2*g * 1j  # *np.sqrt(2)
    return [(-1j * delta - gamma + 2j * alpha * (
            np.abs(B_1) ** 2 + 2 * np.abs(B_2) ** 2)) * B_1 + 2j * alpha * B_2 ** 2 * np.conj(B_1),
            (-1j * delta - gamma + 2 * i_Sigma2_minus_sigma + 2j * alpha * (
                    np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2)) * B_2 +
            2j * alpha * B_1 ** 2 * np.conj(B_2) + p]
# def circuit_1(delta, chi_2, chi_3, gamma_U, gamma_V, gamma_C, sigma, dist, coef, phi, g):
#     def model(z, t, Pi):
#         U_1 = z[0]
#         V_1 = z[1]
#         U_2 = z[2]
#         V_2 = z[3]
#         U_C = z[4]
#         V_C = -gamma_C * U_C + sigma * (U_1 - U_2) + Pi
#         dU_1dt = V_1 - gamma_U * U_1 + (chi_2 * U_1 + chi_3 * U_1 ** 2) * V_1 - V_C
#         dV_1dt = -U_1 - gamma_V * V_1
#         dU_2dt = V_2 - gamma_U * U_2 + (chi_2 * U_2 + chi_3 * U_2 ** 2) * V_2 + V_C
#         dV_2dt = -U_2 - gamma_V * V_2
#         dU_Cdt = -gamma_C * U_C + sigma * (U_1 - U_2) + Pi
#         dzdt = [dU_1dt, dV_1dt, dU_2dt, dV_2dt, dU_Cdt]
#         return dzdt
#
#     # initial condition
#     z_0 = [0.000001, 0, 0.000001, 0, 0]
#
#     # number of time points
#
#     # time points
#     # dist = 10000
#     # coef = 5
#     n = dist * coef  # было умножить
#     t = np.linspace(0, dist, n)
#
#     # step input
#     Pi =2* g * np.cos((1 + delta) * t + phi)  # *np.heaviside(50000-t,0)
#     # (np.exp((1j+delta)*t)).real
#     # store solution
#     U_1 = np.empty_like(t)
#     V_1 = np.empty_like(t)
#     U_2 = np.empty_like(t)
#     V_2 = np.empty_like(t)
#     U_C = np.empty_like(t)
#     # record initial conditions
#     U_1[0] = z_0[0]
#     V_1[0] = z_0[1]
#     U_2[0] = z_0[2]
#     V_2[0] = z_0[3]
#     U_C[0] = z_0[4]
#     # solve ODE
#     for i in range(1, n):
#         # span for next time step
#         tspan = [t[i - 1], t[i]]
#         # solve for next step
#         z = odeint(model, z_0, tspan, args=(Pi[i],))
#         # store solution for plotting
#         U_1[i] = z[1][0]
#         V_1[i] = z[1][1]
#         U_2[i] = z[1][2]
#         V_2[i] = z[1][3]
#         U_C[i] = z[1][4]
#         # next initial condition
#         z_0 = z[1]
#     B_1 = (U_1 + U_2) / 2
#     B_2 = (U_1 - U_2) / 2
#     return t, U_1, U_2
#
