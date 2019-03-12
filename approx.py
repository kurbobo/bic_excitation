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
    Pi =2*g * np.cos((1 +delta) * t + phi)#  *np.heaviside(5000-t,0)
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
    return t, U_1, U_2


def oscillator(delta, gamma, alpha, f,dist,coef):
    def model(z, t):
        X = z[0]
        V = z[1]
        dXdt = V
        dVdt =-gamma*V-X-alpha*X**3+f*np.cos((delta+1)*t)
        dzdt = [dXdt, dVdt]
        return dzdt

    # initial condition
    z_0 = [0, 0]

    # number of time points
    # time points
    # dist = 10000
    # coef = 5
    n = dist * coef  # было умножить
    t = np.linspace(0, dist, n)

    # step input
    #  *np.heaviside(5000-t,0)
    # (np.exp((1j+delta)*t)).real
    # store solution
    X = np.empty_like(t)
    V = np.empty_like(t)
    # record initial conditions
    X[0] = z_0[0]
    V[0] = z_0[1]
    # solve ODE
    for i in range(1, n):
        # span for next time step
        tspan = [t[i - 1], t[i]]
        # solve for next step
        z = odeint(model, z_0, tspan)
        # store solution for plotting
        X[i] = z[1][0]
        V[i] = z[1][1]
        # next initial condition
        z_0 = z[1]
    return t, X, V
def zfunc_osc(z, t, delta, gamma, alpha, f):
    X = z
    print(type(X))
    return [-gamma/2*X - 1j*delta * X + 1.5j*alpha*np.abs(X)**2*X - f/4]
