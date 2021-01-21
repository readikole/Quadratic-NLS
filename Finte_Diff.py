import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import time 

# we are seeking to find a solution to the poisson equation subject to some boundary conditions and an initial cond.
start_time  = time.time()
N = 100
L = 10
x_start = -L
x_end = L  
h = (x_end-x_start)/N         #the total length of 2L
tau = 1
b = 0               #constant coefficient of the nonlinearity


def init(x):
    '''function of initial conditions'''
    y = 1/((np.exp(x-2) + np.exp(-x+2))*(np.exp(x-2) + np.exp(-x+2))/4)
    #y = np.exp(-x**2)
    return y
t = np.zeros(N)
x = np.linspace(x_start, x_end, N)
#define an N*N matrix
u = np.zeros([N, N], dtype = 'cfloat')

for i in range(1, N-1, 1):
    
    for idx in range(1,N-1, 1):
        #intial conditions 
        u[0, 0] = init(x[0])
        u[0, 1] = init(x[1])
        u[i, 1] = (1 + (2j*tau/h**2) + 4j*tau - (1j*b*tau))*u[i, 0] + 1j*(tau/h**2)*(u[i+1, 0] - u[i-1, 0]) - (6j*tau)*u[i,0]**2 - 1j*b*tau*np.conj(u[i, 0])
        u[N-1, 1] = 0
 
for i in range(2, N-1, 1):
    
    for idx in range(2,N-1, 1):
        #intial conditions 
        u[0, 1] = init(x[2])
        u[i, 2] = (1 + (2j*tau/h**2) + 4j*tau - (1j*b*tau))*u[i, 1] + 1j*(tau/h**2)*(u[i+1, 1] - u[i-1, 1]) - (6j*tau)*u[i,1]**2 - 1j*b*tau*np.conj(u[i, 1])
        u[N-1, 2] = 0


for i in range(3, N-1, 1):
    
    for idx in range(2,N-1, 1):
        #intial conditions 
        u[i, idx] = init(x[idx])
        u[i, idx+1] = (1 + (2j*tau/h**2) + 4j*tau - (1j*b*tau))*u[i, idx] + 1j*(tau/h**2)*(u[i+1, idx] - u[i-1, idx]) - (6j*tau)*u[i,idx]**2 - 1j*b*tau*np.conj(u[i, idx])
        u[N-1, idx] = 0
for i in range(N-1):
    t[i] = i*tau

end_time = time.time() - start_time
print('CPU Time: ', end_time)
X, Y = np.meshgrid(x, t)

fig = plt.figure(figsize=(8, 5))
ax = plt.axes(projection='3d')

ax.plot_surface(X, Y, np.real(u),cmap='viridis', edgecolor='none')
ax.plot_surface(X, Y, np.imag(u),cmap='viridis', edgecolor='none')
ax.set_title('Surface plot')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('u')
plt.show()







#3d plotting
fig = plt.figure(figsize=(10, 5))
ax = fig.gca(projection='3d')


'''X, Y = np.meshgrid(x, t)


surf = ax.plot_surface(X, Y, np.real(u), cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()'''
# Customize the z axis
#ax.set_zlim(-1.01, 1.01)
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
#fig.colorbar(surf, shrink=0.5, aspect=5)
