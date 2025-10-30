import numpy as np
import matplotlib.pyplot as plt


## Parameters
pi = 3.1415
N = (0.05 * 1.5)/(pi * (3.5*10**-9)**2)
r = 3.5 * 10**-9 #m
#delP = 15 * 133.3 #mmHg -> Pa (kg/m.s)
delX = 50 * 10**-9 #m

###### Task 2a. ######
## Run Simulation
#Ht = np.linspace(0.1,0.5,100)
#eta = 0.0015 * (1 + 2.5 * Ht) # Pa.s (kg/m.s^2)

###### Task 2b. ######
#Ht = 0.25
#eta = 0.0015 * (1 + 2.5 * Ht) # Pa.s (kg/m.s^2)

#BP = (np.linspace(50, 200, 100))
#delP = (BP - 45) * 133.3 

###### Task 2c. ######
# Create x (BP) and y (Hematocrit) independant variables for z output (filter rate)
Ht = np.linspace(0, 1, 100)
BP = np.linspace(50, 200, 100)
Ht, BP = np.meshgrid(Ht, BP)

# Calculatations
delP = (BP - 45) * 133.3    # Pascals
eta = 0.0015 * (1 + 2.5 * Ht) #Pa.s (kg/m.s^2)

## Simulation ##
Qv = N*(pi * r**4) / (8 * eta) * (delP/delX)

## Show result
plt.figure()
fig = plt.axes(projection = '3d')
fig.contour3D(Ht ,BP ,Qv * 60 * 10**6,100)
fig.set_xlabel('Hematocrit (%)')
fig.set_ylabel('Blood Pressure (mmHg)')
fig.set_zlabel('Microfiltration Rate (mL/min)')
plt.show()



