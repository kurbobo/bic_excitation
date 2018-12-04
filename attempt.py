delta=0.5
sigma=0.102
sigma2=0
gamma=0.06
alpha=0.5
p=0.185
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
if __name__ == "__main__":
    # Define the right-hand-side of the differential equation.
    def zfunc(z, t,gamma,sigma,sigma2,alpha,p):
        B_1,B_2 = z
        return [(-1j*delta-gamma + 2j*alpha*( np.abs(B_1)**2 +2*np.abs(B_2)**2))*B_1 +2j*alpha*B_2**2*np.conj(B_1),
                (-1j * delta - gamma - 2 * sigma + 2j * alpha * (np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2))*B_2 +
        2j * alpha * B_1 ** 2 * np.conj(B_2) + p]
    # Set up the inputs and call odeintz to solve the system.
    z0 = np.array([0+0j, 0.02+0.03j])
    t = np.arange(0, 3000,0.01)
    z, infodict = odeintz(zfunc, z0, t, args=(gamma,sigma,sigma2,alpha,p), full_output=True)
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(t,np.sqrt( z[:,0].real**2+z[:,0].imag**2), label='B_1')
    plt.plot(t,np.sqrt(  z[:,1].real**2+z[:,1].imag**2), label='B_2')
    plt.xlabel('t')
    plt.grid(True)
    plt.legend(loc='best')
    plt.show()
