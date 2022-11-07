# Homework 4 🥵🥵
B08504055 田沛茵🈹

* Necessary imports and parameters
```python {.line-numbers}
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
```

* Analytical solution to the diffusion equation (n=50)
```python
def C(x,t):
    a, b = 0, 0
    for n in range(1,51):
        a += np.cos(n*pi*x/L)/n**2*np.exp(-D*t*n**2*pi**2/L**2)
        b += (-1)**n*np.cos(n*pi*x/L)/n**2*np.exp(-D*t*n**2*pi**2/L**2)
        
    return te*k2*(D*t/L**2+(L-x)**2/2/L**2-1/6-2*a/pi**2) + \
            (1-te)*k2*(D*t/L**2+x**2/2/L**2-1/6-2*b/pi**2)
```

### (a) Concentration profile with $t_e = 0.001$
![](https://i.imgur.com/GL6bkx2.png)

```python
te = 0.001
for i in range(7):
    plt.plot(x, C(x,t[i]), label='t = '+str(t[i])+'s')
```

### (b) U vs. t with $t_e = 0.001$
![](https://i.imgur.com/ECaHnvk.png)

```python
tx = np.linspace(0,3600,1000)
sigma = 96500**2*D/(k1*te*(1-te))
U = -I*L/sigma-te/96500*k1*C(0,tx)-(1-te)/96500*k1*C(L,tx)

plt.plot(tx,U)
```

### (c) U vs. $\sqrt{t}$ with $t_e = 0.001$
![](https://i.imgur.com/yabW5zR.png)

```python
t_sqrt = np.sqrt(tx)

# proportional line
slope = (U[100]-U[0])/(t_sqrt[100])
U_approx = slope*t_sqrt + U[0]

plt.plot(t_sqrt,U,label='U')
plt.plot(t_sqrt,U_approx,label='Approximation',linestyle='dashed',color='black')
```

$t_1 \approx 625$

$\tau^\delta=\dfrac{L^2}{D^\delta}=2500$

$\Rightarrow t_1 \approx \dfrac{1}{4} \tau^\delta$

### (d) Concentration profile with $t_e = 0.999$
![](https://i.imgur.com/Jtk4mIG.png)

```python
te = 0.999
for i in range(7):
    plt.plot(x, C(x,t[i]), label='t = '+str(t[i])+'s')
```

### (e) U vs. t with $t_e = 0.999$
![](https://i.imgur.com/2lEYdxX.png)

```python
sigma = 96500**2*D/(k1*te*(1-te))
U = -I*L/sigma-te/96500*k1*C(0,tx)-(1-te)/96500*k1*C(L,tx)

plt.plot(tx,U,label='te=0.999',linestyle='dashed')
```

The two curves are identical
