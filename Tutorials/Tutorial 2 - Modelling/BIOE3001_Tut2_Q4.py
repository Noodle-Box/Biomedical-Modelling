import numpy as np
import matplotlib.pyplot as plt


## Time course
dt = 1 #s
T = 60*60*10 #s

## Parameters
pi = 3.1415
r = 2*10**-10 #m
Npores = (0.05 * 1.5)/(pi * (r)**2)
eta = 0.002 # Pa.s (kg/m.s^2)
delX = 50*10**-9 #m
Vbladder = 20; #mL per kidney
KLN = 0.00028 #1/s
rblood = 0.0016 #mg/mL.s
PeePct = 0.9 #frac of full

## Initial Conditions
V = [None] * (T//dt)
N = [None] * (T//dt)
C = [None] * (T//dt)
V[0] = 1.0150 #mL
N[0] = 4 #mg
C[0] = 8 #mg/mL

## Setting Pee Rate
PeeFreq = np.linspace(0,T//dt,10) #s

## Run Simulation
for t in range(0,(T//dt-1)):
    if t < (T//dt)*0.5:
        delP = (60-45)*133.3 
    else:
        delP = (120-45)*133.3
    PeeIn = Npores*(pi*r**4)/(8*eta)*(delP/delX)*10**6
    if sum(np.equal(t,PeeFreq)):
        dV = PeeIn*(1-V[t]/Vbladder) - V[t]*PeePct
        dN = KLN*(C[t]-N[t]/V[t])*V[t] - N[t]*PeePct
    else:
        dV = PeeIn*(1-V[t]/Vbladder)
        dN = KLN*(C[t]-N[t]/V[t])* V[t]
    dC = rblood - KLN*(C[t]-N[t]/V[t])
    V[t+1] = V[t] + dV
    N[t+1] = N[t] + dN
    C[t+1] = C[t] + dC

## Plot Figure
plt.figure()
plt.subplot(311)
plt.plot(range(0,T//dt), C)
plt.ylabel('Blood Urea (g/mL)')
plt.subplot(312)
plt.plot(range(0,T//dt), V)
plt.ylabel('Urine (mL)')
plt.subplot(313)
plt.plot(range(0,T//dt), N)
plt.xlabel('Time (s)')
plt.ylabel('Urea (g)')
plt.show()


