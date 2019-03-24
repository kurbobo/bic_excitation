import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
fig = plt.figure()
ax = plt.axes(projection='3d')
delta=0.2#0.2581988897471611
alpha=1
delta_light=0.2
gamma_light=0.4
gamma=0.1#0.11547005383792516
y=np.arange(-1000,1000,0.001)
B_2_abs=np.sqrt(gamma**2/4/alpha**2/y**2)
b_abs=np.sqrt(delta/alpha-3*gamma**2/4/alpha**2/y**2-y**2)
p_1=(delta-delta_light)*y-gamma**2/alpha/y
p_2=gamma**3/alpha**2/y**3+2*gamma*y -gamma/2/alpha/y*(3*delta-delta_light)
p_abs=np.sqrt(p_1**2+p_2**2)
ax.plot3D(b_abs,B_2_abs,p_abs, 'blue')
plt.show()