delta=0
sigma=0.0317010832
sigma2=0.1
gamma=0.06
alpha=0
p=0.185
dist=20000
FD=26000
import numpy as np
from numpy.fft import *
from scipy.integrate import odeint
import real_curcuit as rc
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
                (-1j * delta - gamma - 2 * sigma +2j * sigma2 + 2j * alpha * (np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2))*B_2 +
        2j * alpha * B_1 ** 2 * np.conj(B_2) + p]
    # Set up the inputs and call odeintz to solve the system.
    z0 = np.array([0.00001, 0])
    t = np.arange(0, dist,1)
    # real_curcuit(delta, chi_2, chi_3, R, r,zed, R_, C,  varepsilon_0, dist, coef):
    t_1,B_1r,B_2r=rc.real_curcuit(delta,#delta
                                  0,#chi2
                                  0,#chi_3
                                  0.10897247358851686,#R
                                  0.1,#r
                                  0.1,#zed
                                  0.002,#R_
                                  4,#k=C_/C
                                  0.142552073818657,#varepsilon
                                  dist,#dist
                                  1,#coef
                                  fi)#phi
    for fi in np.arange(0, 2 * 3.14, 0.2):
        z, infodict = odeintz(zfunc, z0, t, args=(gamma,sigma,sigma2,alpha,p), full_output=True)
    # import matplotlib.pyplot as plt
    # plt.clf()
    # #print slow amplitudes
    # # plt.plot(t,np.abs(z[:,0]), label='B_1')
    # plt.plot(t,np.abs(z[:,1]), label='B_2')
    # #print real values
    # # plt.plot(t_1, np.abs(B_1r), label='B_1 real')
    # # plt.plot(t_1, np.abs(B_2r), label='B_2 real')
    # plt.xlabel('t')
    # plt.grid(True)
    # plt.legend(loc='best')
    # plt.show()
