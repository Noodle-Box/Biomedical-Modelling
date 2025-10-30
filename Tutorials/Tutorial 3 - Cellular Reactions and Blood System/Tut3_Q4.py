import numpy as np
import matplotlib.pyplot as plt

## Parameters
k1 = 1      # 1/uM*s
k2 = 0.15   #1/s
CEO = 5     #uM

Kmax = CEO * k2


## Independent Variables
CS, k_1 =   np.meshgrid(np.linspace(0,40,100),np.linspace(0,40,100)) #uM, 1/s

Km = (k_1 + k2)/k1

                
## Dependent Equation and Plot
r = Kmax * CS/(Km + CS)

ax = plt.axes(projection='3d')
ax.contour3D(CS, k_1, r, 100)
ax.set_xlabel('Substrate, uM')
ax.set_ylabel('k1neg, uM')
ax.set_zlabel('Rate, uM/s');
ax.set_title("MM CS vs k_1 simultaneous change")
plt.show()
