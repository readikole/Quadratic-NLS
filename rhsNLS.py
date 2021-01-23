import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import axes3d
import matplotlib.cm as cm

#intialize constants
b = 0            
L = 60
N = 1000
dx = L/N
x = np.arange(-L/2, L/2, dx)

#Define discrete wave numbers
kappa = 2*np.pi*np.fft.fftfreq(N, d=dx)

#initial condition
u0 = 1/np.cosh(x)**2


#simulate PDE
dt = 0.01
t = np.arange(0, 55500*dt, dt)


#function to be integrated
def rhsNLS(u,t , kappa, b):
    uhat = np.fft.fft(u)
    #d_uhat = (1j)*kappa*uhat
    dd_uhat =  - np.power(kappa, 2)*uhat
    #d_u  = np.fft.ifft(d_uhat)
    dd_u = np.fft.ifft(dd_uhat)
    du_dt = - (1j)*(dd_u + 6*u*u - 4*u +b*(u - np.conj(u)))
    return du_dt.real

u = odeint(rhsNLS, u0, t, args=(kappa, b))




#waterfall plot
fig  =  plt.figure()
ax = fig.add_subplot(111, projection='3d')
u_plot = u[0:-1:10, :]
plt.set_cmap('jet_r')
for j in range(u_plot.shape[0]):
    ys = j*np.ones(u_plot.shape[1])
    ax.plot(x, ys, u_plot[j, :], color = cm.jet(j*20))
plt.show()


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
