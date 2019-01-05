from sympy import *
# B_1=Symbol('B_1',complex=True,)
# B_2=Symbol('B_2',complex=False)
from matplotlib import pyplot as plt
import numpy as np
delta=0.5
sigma=0.102
gamma=0.06
alpha=0.5
p=0.185
# gamma,alpha=symbols('gamma alpha',complex=False)
# answ=(simplify(solve(-gamma+2j*alpha*(B_1*conjugate(B_1)+2*B_2**2)+2j*alpha*B_2**2*conjugate(B_1)/B_1,B_1)))
# answer=0.5*np.sqrt(-12.0*B_2**2 - 2.0*I*gamma/alpha)
def fun(B_2):
    answ=0.5*np.sqrt(-12.0*B_2**2 - 2.0*1j*gamma/alpha)
    return np.abs(answ)
# B_2=np.arange(0.01,p,0.01)
# plt.plot(B_2,fun(B_2))
# plt.show()
print(fun(0.41))

