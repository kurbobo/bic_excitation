from matplotlib import pyplot as plt
import numpy as np
Gamma=0.01
sigma=0.01
sigma2=0.1
delta_0=0.1
alpha=0.01
B=np.arange(0.1,5,0.1)
p=B*(Gamma+2*sigma)+1j*(2*sigma2+delta_0-2*alpha*B**2)*B
# print()
plt.plot(B,np.abs(p))
xplus=np.sqrt((2*(2*sigma2+delta_0)+np.sqrt((2*sigma2+delta_0)**2-3*(Gamma+2*sigma)**2))/(6*alpha))
xminus=np.sqrt((2*(2*sigma2+delta_0)-np.sqrt((2*sigma2+delta_0)**2-3*(Gamma+2*sigma)**2))/(6*alpha))
B_minus=np.sqrt((4*delta_0-2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)
B_plus=np.sqrt((4*delta_0+2*np.sqrt(4*delta_0**2-3*Gamma**2))/alpha)
plt.axvline(x=xplus,color='red')
plt.axvline(x=xminus,color='red')
plt.axvline(x=B_plus,color='green')
plt.axvline(x=B_minus,color='green')
plt.ylabel('|p|')
plt.xlabel('B_{20}')
plt.text(xminus+0.1,5,'$x_{-}$',fontsize=14,color='red')
plt.text(xplus+0.1,5,'$x_{+}$',fontsize=14,color='red')
plt.text(B_minus+0.1,5,'$y_{-}$',fontsize=14,color='green')
# plt.text(B_plus+0.1,5,'$y_{+}$',fontsize=14,color='green')
plt.show()
