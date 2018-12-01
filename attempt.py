delta=1
gamma=2
sigma=0.8
sigma2=1
alpha=0.01
p=1
C=-delta-2*sigma2
A=-gamma-2*sigma
# x=(A* p)/((A**2+C**2)**3)* (-(A**2+C**2)**2+4*alpha *C*p**2)
# y=(C* p* (A**2+C**2)**2+2*alpha*p**3*(A**2-C**2))/(A**2+C**2)**3
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
        return [-(gamma+1j*delta)*B_1+2j*alpha*(np.abs(B_1)**2+2*np.abs(B_2)**2)*B_1 + 2j*alpha*B_2**2*np.conj(B_1),
                -(gamma+1j*delta) * B_2-2*(sigma-1j*sigma2)*B_2 + 2j * alpha * (
                            np.abs(B_2) ** 2 + 2 * np.abs(B_1) ** 2) * B_2 + 2j * alpha * B_1 ** 2 * np.conj(B_2)+p]

    # Set up the inputs and call odeintz to solve the system.
    z0 = np.array([0+0j, 0+0j])
    t = np.linspace(0, 10, 101)
    z, infodict = odeintz(zfunc, z0, t, args=(gamma,sigma,sigma2,alpha,p), full_output=True)
    import matplotlib.pyplot as plt
    plt.clf()
    plt.plot(t,np.sqrt( z[:,0].real**2+z[:,0].imag**2), label='B_1')
    plt.plot(t,np.sqrt(  z[:,1].real**2+z[:,1].imag**2), label='B_2')
    a=((np.sqrt(z[:,1].real**2+z[:,1].imag**2))[100])

    #------------------------------------------------------------
    x_0 = -A * p / (A ** 2 + C ** 2)
    y_0 = C * p / (A ** 2 + C ** 2)
    x_1 = p ** 3 * (A ** 2 - C ** 2) / (A ** 2 + C ** 2) ** 3
    y_1 = 2 * p ** 3 * (A ** 2 - C ** 2) / (A ** 6 + 3 * A ** 4 * C ** 2 + 3 * A ** 2 * C ** 4 + C ** 6)
    x_2 = 2 * p ** 5 * (
                2 * A ** 5 - 5 * A ** 4 * C + 8 * A ** 3 * C ** 2 + 4 * A ** 2 * C ** 3 - 10 * A * C ** 4 + C ** 5) / (
                      A ** 12 + 6 * A ** 10 * C ** 2 + 15 * A ** 8 * C ** 4 + 20 * A ** 6 * C ** 6 + 15 * A ** 4 * C ** 8 + 6 * A ** 2 * C ** 10 + C ** 12)
    y_2 = p ** 5 * (
                -6 * A ** 5 + 4 * A ** 4 * C + 8 * A ** 3 * C ** 2 - 16 * A ** 2 * C ** 3 - 2 * A * C ** 4 + 12 * C ** 5) / (
                      A ** 12 + 6 * A ** 10 * C ** 2 + 15 * A ** 8 * C ** 4 + 20 * A ** 6 * C ** 6 + 15 * A ** 4 * C ** 8 + 6 * A ** 2 * C ** 10 + C ** 12)
    x = x_0
    y = y_0
    print('B_20=', np.sqrt(x ** 2 + y ** 2), '0-й порядок')
    x = x_0 + alpha * x_1
    y = y_0 + alpha * y_1
    print('B_20=', np.sqrt(x ** 2 + y ** 2), '1-й порядок')
    x = x_0 + alpha * x_1 + alpha ** 2 * x_2
    y = y_0 + alpha * y_1 + alpha ** 2 * y_2
    print('B_20=', np.sqrt(x ** 2 + y ** 2), '2-й порядок')
    print(a,'- то что должно быть')
    plt.xlabel('t')
    plt.grid(True)
    plt.legend(loc='best')
    plt.show()
