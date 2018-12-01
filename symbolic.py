import numpy as np
from matplotlib.pyplot import *
from sympy import *
from sympy.solvers import solve
A=Symbol('A')
p=Symbol('p')
C=Symbol('C')
alpha=Symbol('alpha')
x_2,y_2= symbols('x_2 y_2')
x_0=-A*p/(A**2+C**2)
y_0=C*p/(A**2+C**2)
x_1=p**3*(A**2 - C**2)/(A**2 + C**2)**3
y_1=2*p**3*(A**2 - C**2)/(A**6 + 3*A**4*C**2 + 3*A**2*C**4 + C**6)
x_2=2*p**5*(2*A**5 - 5*A**4*C + 8*A**3*C**2 + 4*A**2*C**3 - 10*A*C**4 + C**5)/(A**12 + 6*A**10*C**2 + 15*A**8*C**4 + 20*A**6*C**6 + 15*A**4*C**8 + 6*A**2*C**10 + C**12)
y_2=p**5*(-6*A**5 + 4*A**4*C + 8*A**3*C**2 - 16*A**2*C**3 - 2*A*C**4 + 12*C**5)/(A**12 + 6*A**10*C**2 + 15*A**8*C**4 + 20*A**6*C**6 + 15*A**4*C**8 + 6*A**2*C**10 + C**12)
x=x_0+alpha*x_1
y=y_0+alpha*y_1
print('B_20=',(x**2+y**2),'1-й порядок')
x=x_0+alpha*x_1+alpha**2*x_2
y=y_0+alpha*y_1+alpha**2*y_2
print('B_20=',(x**2+y**2),'2-й порядок')
# eqs=(-A*x_1+C*y_1+2*x_0**2*y_0+2*y_0**3,
# A*y_1+C*x_1+2*x_0*y_0**2+2*x_0**3)

# eqs=(-A*x_2+C*y_2+4*y_0*(x_0*x_1+y_0*y_1)+2*y_1*(x_0**2+y_0**2),
#             A * y_2 + C * x_2 + 4 * x_0 * (x_0 * x_1 + y_0 * y_1) + 2 * x_1 * (x_0 ** 2 + y_0 ** 2))
# print(simplify(solve(eqs,x_2,y_2)))
