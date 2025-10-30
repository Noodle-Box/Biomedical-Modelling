import numpy as np
import matplotlib.pyplot as plt

## Parameters 

F = 9.649*10**4     # (C/mol)
R = 8.314           # (J/mol*k)
T = 310             # (Kelvin, Temp)
z = 2

# Concentrations of Sodium (Ca2+)
Co_Ca = 1.2*10**-3   # (M)
Ci_Ca = 0.1*10**-6    # (M)

E_Ca = ((R * T) / (z * F)) * np.log(Co_Ca / Ci_Ca)

print(f"Calcium:      E = {E_Ca:.4f} V\n")
