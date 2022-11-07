# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 00:48:51 2022

@author: sprit
"""

# imports
import numpy as np
import matplotlib.pyplot as plt

# parameters
I = 0.01 # A/cm^2
D = 1e-4 # cm^2/s
L = 0.5 # cm
k1 = 1e6 # J/mol
k2 = I*L/96500/D
t = [60, 300, 600, 1200, 1800, 2700, 3600]
x = np.linspace(0,0.5,100)
pi = np.pi

def C(x,t):
    a, b = 0, 0
    for n in range(1,51):
        a += np.cos(n*pi*x/L)/n**2*np.exp(-D*t*n**2*pi**2/L**2)
        b += (-1)**n*np.cos(n*pi*x/L)/n**2*np.exp(-D*t*n**2*pi**2/L**2)
        
    return te*k2*(D*t/L**2+(L-x)**2/2/L**2-1/6-2*a/pi**2) + \
            (1-te)*k2*(D*t/L**2+x**2/2/L**2-1/6-2*b/pi**2)
            

# (a)
te = 0.001

plt.figure(0)
for i in range(7):
    plt.plot(x, C(x,t[i]), label='t = '+str(t[i])+'s')
    
plt.xlabel('L')
plt.ylabel('Concentration')
plt.legend(fancybox=True, framealpha=0.5, fontsize=8)
plt.savefig('hw4_a.png', dpi=200)


# (b)
tx = np.linspace(0,3600,1000)
sigma = 96500**2*D/(k1*te*(1-te))
U = -I*L/sigma-te/96500*k1*C(0,tx)-(1-te)/96500*k1*C(L,tx)

plt.figure(1)
plt.plot(tx,U)
plt.xlabel('t')
plt.ylabel('U')
plt.savefig('hw4_b.png', dpi=200)

plt.figure(4)
plt.plot(tx,U,label='te = 0.001')
plt.xlabel('t')
plt.ylabel('U')


# (c)
t_sqrt = np.sqrt(tx)

# proportional line
slope = (U[100]-U[0])/(t_sqrt[100])
U_approx = slope*t_sqrt + U[0]

plt.figure(2)
plt.plot(t_sqrt,U,label='U')
plt.plot(t_sqrt,U_approx,label='Approximation',linestyle='dashed',color='black')
plt.plot(25, slope*25+U[0], 'rx')
plt.text(27, slope*25+U[0], r'$\sqrt{t_1} \approx 25$', color='red')
plt.xlabel(r'$\sqrt{t}$')
plt.ylabel('U')
plt.legend()
plt.savefig('hw4_c.png', dpi=200)


# (d)
te = 0.999

plt.figure(3)
for i in range(7):
    plt.plot(x, C(x,t[i]), label='t = '+str(t[i])+'s')
    
plt.xlabel('L')
plt.ylabel('Concentration')
plt.legend(fancybox=True, framealpha=0.5, fontsize=8)
plt.savefig('hw4_d.png', dpi=200)


# (e)
sigma = 96500**2*D/(k1*te*(1-te))
U = -I*L/sigma-te/96500*k1*C(0,tx)-(1-te)/96500*k1*C(L,tx)

plt.figure(4)
plt.plot(tx,U,label='te = 0.999',linestyle='dashed')
plt.xlabel('t')
plt.ylabel('U')
plt.legend(fancybox=True, framealpha=0.5, fontsize=8)
plt.savefig('hw4_e.png', dpi=200)

