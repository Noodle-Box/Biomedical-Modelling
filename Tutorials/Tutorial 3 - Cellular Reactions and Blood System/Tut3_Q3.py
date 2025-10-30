import numpy as np
import matplotlib.pyplot as plt

"MM Kinetics in week 1 urine Peeeee example"

## Time course
dt = 1 #s
T = 60*60*10 #s

## Parameters
Vbladder = 20; #mL per kidney
PeeIn = 0.015 #mL/s
KLN = 0.00028 #mg/s
PeePct = 0.9 #frac of full

## MM Parameters
k_1 = 0.003     # mL/mg*s
k1 = 0.001      # 1/s
k2 = 0.03         # 1/s
CE0 = 0.1       # mg/mL

Kmax = CE0 * k2
Km = (k_1 + k2)/k1

## Initial Conditions
V = [None] * (T//dt)
N = [None] * (T//dt)
C = [None] * (T//dt)

V[0] = 1.0150 #mL
N[0] = 4 #mg
C[0] = 8 #mg/mL #changed in wk2

## Setting Pee Rate
PeeFreq = np.linspace(0,T//dt,10) #s

## Run Simulation
for t in range(0,(T//dt-1)):
    if sum(np.equal(t,PeeFreq)):
        dV = PeeIn*(1-V[t]/Vbladder) - V[t]*PeePct
        dN = KLN*(C[t]-N[t]/V[t])* Vbladder - N[t]*PeePct
    else:
        dV = PeeIn*(1-V[t]/Vbladder)
        dN = KLN*(C[t]-N[t]/V[t])*Vbladder
    dC = ((Kmax * C[t])/(Km + C[t])) - KLN*(C[t]-N[t]/V[t])
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
