# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 12:01:36 2022

@author: u92780tw
"""

import numpy as np
from scipy.fft import fft, ifft, fftshift
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

dx = 0.01 #space step
L = 100 #length of x domain
c = -4.0 #non-linear coefficient
iterations = 2000
divisor = 50
dt = 0.01 #time step

def IC(x):
    return 2 / np.cosh(2*(x)) #+ 0.8 / np.cosh(0.8*(x-1))

def split_step(dx, L, dt, c, iterations, IC, divisor):

    N = int(L/dx)  #number of discrete space steps
    x = np.arange(-L/2, L/2, L/N)
    k = (2*np.pi/L)*fftshift(np.arange(-N/2,N/2)) # k_space
    y = IC(x)
    y_array = y ## initializes an array to store y values for each time step
    final_time = iterations*dt

    for i in range(1, iterations+1):

        #EXECUTES THE FIRST POTENTIAL STEP
        y = np.exp(-0.5 * c * np.abs(y)**2 * dt * 1j, dtype=complex) * y

        #FOURIER TRANSFORMS TO MOMENTUM SPACE
        y = fft(y)

        #EXECUTES THE FIRST KINETIC STEP
        y = np.exp(-(k ** 2) * dt * 1j, dtype=complex) * y

        #INVERSE FOURIER TRANSFORMS FROM MOMENTUM SPACE
        y = ifft(y)

        #EXECUTES THE SECOND POTENTIAL STEP
        y = np.exp(-0.5 * c * np.abs(y)**2 * dt * 1j, dtype=complex) * y

        if i%(iterations/divisor) == 0:
            y_array = np.column_stack((y_array, np.abs(y)))

    return x, y_array

def RK4(dx, L, dt, lamda, iterations, IC, divisor):

    N = int(L/dx)  #number of discrete space steps
    x = np.arange(-L/2, L/2, L/N)
    k = (2*np.pi/L)*fftshift(np.arange(-N/2,N/2))
    sigma = -0.5
    ik2sigma = 1j* k**2 * sigma

    u = IC(x)
    udata = np.abs(u)
    U = fft(u)
    for i in range(1, iterations+1):
        t = i*dt
        g = 1j*lamda*dt
        Ep = np.exp(dt*ik2sigma/2)
        Ep2 = Ep**2;
        Em = np.exp(-dt*ik2sigma/2)
        Em2 = Em**2;
        Etp = np.exp(t*ik2sigma);
        Etm = np.exp(-t*ik2sigma)

        a = g*Etp*fft(np.abs(ifft(Etm*U))**2 * (ifft(Etm*U)))
        b = g*Etp*Ep*fft(np.abs(ifft(Etm*Em*(U+a/2)))**2 * (ifft(Etm*Em*(U+a/2))))
        c = g*Etp*Ep*fft(np.abs(ifft(Etm*Em*(U+b/2)))**2 * (ifft(Etm*Em*(U+b/2))))
        d = g*Etp*Ep2*fft(np.abs(ifft(Etm*Em2*(U+c)))**2 * (ifft(Etm*Em2*(U+c))))
        U = U + (a + 2*(b+c) + d)/6

        if i%(iterations/divisor) == 0:
            udata = np.column_stack((udata, np.abs(ifft(U*Etm))))

    return x, udata

x, y_array = split_step(dx, L, dt, c, iterations, IC, divisor)


y_array = np.delete(y_array, slice(0,4600), axis=0)
y_array = np.delete(y_array, slice(800, 5400), axis=0)

x = np.delete(x, slice(0,4600), axis=0)
x = np.delete(x, slice(800, 5400), axis=0)

t = np.arange(0, y_array.shape[1], 1) * dt * iterations / 50
X, T = np.meshgrid(x, t)

Z = y_array.transpose()

fig = plt.figure()
ax = plt.axes(projection='3d')

ax.plot_surface(X, T, Z, rstride=2, cstride=2, cmap='viridis_r')
#ax.contour3D(X, Y, Z, 100, cmap='binary')

print(y_array[400][0]-y_array[400][50])

ax.view_init(25, 120)

ax.set_xlabel('z')
ax.set_ylabel('time')
ax.set_zlabel('U(z,t)')

plt.tight_layout(pad=0, w_pad=0, h_pad=0)
plt.savefig("figure.eps", dpi=1500, format="eps")

plt.show()


