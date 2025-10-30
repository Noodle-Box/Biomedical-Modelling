import numpy as np
import matplotlib.pyplot as plt


"MM Kinetics of Concentration over time"

# Time paramters

T = 100 #s
dt = 1 #s
CS = [None] * (T//dt)
CS[0] = 10

# Equation parameters (SI units)

Kmax = 5*0.15
Km = (20 + 0.15)/5

# Simulate
for t in range(0, T//dt-1):
    
    # Compute MM
    R = (Kmax * CS[t])/(Km + CS[t])
    
    # Add to "next step"
    CS[t+1] = CS[t] - R
    
## Dependent Equation and Plot
plt.figure()
plt.plot(range(0,(T//dt)),CS,color='k')
plt.title('MM Kinetics')
plt.ylabel('Concentration (uM)')
plt.xlabel('Time (s)')



    
