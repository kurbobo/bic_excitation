import numpy as np
delta=2
gamma=1
sigma=0.8
sigma2=1
alpha=0.01
p=2
C=-delta-2*sigma2
A=-gamma-2*sigma
x = -A * p / (A ** 2 + C ** 2)
y = C * p / (A ** 2 + C ** 2)
B_20=np.sqrt(x**2+y**2)
print('B_20=',B_20)
if(delta<=2*gamma and delta>=np.sqrt(3)*gamma):
    print('delta<=2*gamma and delta>=np.sqrt(3)*gamma')
    print('sqrt((delta-2*np.sqrt(delta**2-3*gamma**2))/alpha)=',
          np.sqrt((delta-2*np.sqrt(delta**2-3*gamma**2))/alpha))
    print('sqrt((delta+2*np.sqrt(delta**2-3*gamma**2))/alpha)=',
          np.sqrt((delta + 2 * np.sqrt(delta ** 2 - 3 * gamma ** 2)) / alpha))
    if (B_20>np.sqrt((delta-2*np.sqrt(delta**2-3*gamma**2))/alpha) and B_20<np.sqrt((delta+2*np.sqrt(delta**2-3*gamma**2))/alpha)):
        print('bic is excited!')
    else: print("bic isn't excited")
elif(delta>2*gamma):
    print('delta>2*gamma')
    print('sqrt((delta+2*np.sqrt(delta**2-3*gamma**2))/alpha)=',
          np.sqrt((delta + 2 * np.sqrt(delta ** 2 - 3 * gamma ** 2)) / alpha))
    if(B_20>0 and B_20<np.sqrt((delta+2*np.sqrt(delta**2-3*gamma**2))/alpha)):
        print('bic is excited!')
    else: print("bic isn't excited")
else:
    print("fund bic isn't excited")


