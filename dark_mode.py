import numpy as np
from matplotlib.pyplot import *
from sympy import *
from sympy.solvers import solve
x, A,B,p,alpha = symbols("x A B p alpha")
print(solve(x*sqrt(A**2+(B-2*alpha*x**2)**2)-p,x))
