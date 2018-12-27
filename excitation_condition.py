import numpy as np
delta=0.01
gamma= 0.0018973665961010276
alpha= 0.010349999999999998
p= 0.12023303494693738
sigma= 0.00730220604733538
sigma2= 0.070992366412213748
B=delta+2*sigma2
A=gamma+2*sigma
coef=[4*alpha**2,2*alpha*B,A**2+B**2,-p**2]
roots=np.roots(coef)
print(roots)
for i in roots:
    if(i.imag!=0):
        print('bijection')
        break
# if(delta>np.sqrt(3/4)*gamma):
#     print('delta>np.sqrt(3/4)*gamma')
#     if(B**2-3*A**2<0):
#         print('no extrem, excites')
#     if (B ** 2 - 3 * A ** 2 ==0):
#         print('whata fuck')
#     if (B ** 2 - 3 * A ** 2 > 0):
#         y_minus=np.sqrt((4*delta-2*np.sqrt(4*delta**2-3*gamma**2))/alpha)
#         y_plus = np.sqrt((4 * delta + 2 * np.sqrt(4 * delta ** 2 - 3 * gamma ** 2)) / alpha)
#         print('y_minus=',y_minus,', y_plus=',y_plus)
#         x_minus = np.sqrt(2*B-np.sqrt(B**2-3*A**2) / alpha)
#         x_plus = np.sqrt(2*B+np.sqrt(B**2-3*A**2) / alpha)
#         print('x_minus=', x_minus, ', x_plus=', x_plus)
#         if(y_minus<x_minus and x_plus>y_plus):
#             if()
# else:
#     print("fund bic isn't excited")



