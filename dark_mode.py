import numpy as np
from matplotlib.pyplot import *
from sympy import *
from sympy.solvers import solve
A=Symbol('A')
p=Symbol('p')
C=Symbol('C')
alpha=Symbol('alpha')
B_01,B_02,g,D,x,y= symbols('B_01 B_02 g D x y')
x0=-A*p/(A**2+C**2)
y0=C*p/(A**2+C**2)
B_01=-A*p/(A**2 + C**2) + alpha*p**3*(A**2 - C**2)/(A**2 + C**2)**3
B_02=C*p/(A**2 + C**2) + 2*alpha*p**3*(A**2 - C**2)/(A**6 + 3*A**4*C**2 + 3*A**2*C**4 + C**6)
eqs=(-g*x-D*y-4*alpha*x*B_01*B_02+2*alpha*y*(B_01**2-B_02**2),
-g*y+D*x+2*alpha*y*B_01*B_02+2*alpha*x*(B_01**2-B_02**2))
print(simplify(solve(eqs,x,y)))
