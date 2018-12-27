from matplotlib import pyplot as plt
import numpy as np
Gamma= 6e-04*1e02
alpha= 1e-02*1e02
sigma= 0.12e-02*1e02
sigma2=0# 4.999995000005001e-08
delta_0= 0.5e-02*1e02
B=np.arange(0.01,15,0.001)
p=B*(Gamma+2*sigma)+1j*(2*sigma2+delta_0-2*alpha*B**2)*B
plt.plot(np.abs(p),B)
xplus=np.sqrt((2*(2*sigma2+delta_0)+np.sqrt((2*sigma2+delta_0)**2-3*(Gamma+2*sigma)**2))/(6*alpha))
xminus=np.sqrt((2*(2*sigma2+delta_0)-np.sqrt((2*sigma2+delta_0)**2-3*(Gamma+2*sigma)**2))/(6*alpha))
B_minus=np.sqrt((4*delta_0-2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)
B_plus=np.sqrt((4*delta_0+2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)
plt.axhline(y=xplus,color='red')
plt.axhline(y=xminus,color='red')
plt.axhline(y=B_plus,color='green')
plt.axhline(y=B_minus,color='green')
plt.xlabel('|p|')
plt.ylabel('B_{20}')
# plt.text(xminus+0.1,5,'$x_{-}$',fontsize=14,color='red')
# plt.text(xplus+0.1,5,'$x_{+}$',fontsize=14,color='red')
# plt.text(B_minus+0.1,5,'$y_{-}$',fontsize=14,color='green')
# plt.text(B_plus+0.1,5,'$y_{+}$',fontsize=14,color='green')
plt.show()
