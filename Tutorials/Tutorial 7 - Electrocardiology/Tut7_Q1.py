import numpy as np
import matplotlib.pyplot as plt

## Time course
dt = 1 #ms
T = 300 #s

## Setting Channel Openings
gK = [0]*(T//dt)
gNa = [0]*(T//dt)
gCa = [0]*(T//dt)

Em = [0]*(T//dt)

I_Na = [0]*(T//dt)
I_K = [0]*(T//dt)
I_Ca = [0]*(T//dt)

## Em Parameters
E_Na = 0.071
E_K = -0.094
E_Ca = 0.125

## Run Simulation
for t in range(0,(T//dt)):
    gNa[t] = 5 + 500*(100/np.sqrt(2*3.1415*50**2)*np.exp(-((t-20)**2)/(2*9**2)))
    gCa[t] = 1 + 500*(200/np.sqrt(2*3.1415*130**2)*np.exp(-((t-80)**2)/(2*40**2)))
    gK[t] = 500-10**4.4/(t*1*np.sqrt(2*3.1415))*np.exp(-((np.log(t)-4)**2)/(2*1**2))
    g_sum = gNa[t] + gCa[t] + gK[t]
    
    Em[t] = ((gNa[t]*E_Na)/g_sum) + ((gK[t]*E_K)/g_sum) + ((gCa[t]*E_Ca)/g_sum)
    
    # Ion specific currents
    I_Na[t] = gNa[t]*(Em[t]-E_Na)
    I_K[t] = gK[t]*(Em[t]-E_K)
    I_Ca[t] = gCa[t]*(Em[t]-E_Ca)
    
## Plot Figure
plt.figure()
plt.title("Action Potential")
plt.plot(range(0,T//dt), Em, color='k')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential, (mV)')

plt.figure()
plt.title("Ion specific Conductance")
plt.plot(range(0,T//dt), gNa, color='r')
plt.plot(range(0,T//dt), gCa, color='b')
plt.plot(range(0,T//dt), gK, color='k')
plt.xlabel('Time (ms)')
plt.ylabel('Conductance (mA)')
plt.legend(['gNa','gCa','gK'])

plt.figure()
plt.title("Ion specific current")
plt.plot(range(0,T//dt), I_Na, color='r')
plt.plot(range(0,T//dt), I_Ca, color='b')
plt.plot(range(0,T//dt), I_K, color='k')
plt.xlabel('Time (ms)')
plt.ylabel('Ion Current (mA)')
plt.legend(['gNa','gCa','gK'])
plt.show()

