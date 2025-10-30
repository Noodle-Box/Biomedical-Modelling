import numpy as np
import matplotlib.pyplot as plt

## Parameters 

F = 9.649*10**4     # (C/mol)
R = 8.314           # (J/mol*k)
T = 310             # (Kelvin, Temp)

z1 = 1
z2 = 2

# Concentrations of Sodium (Na+)
Co_Na = 145*10**-3   # (M)
Ci_Na = 10*10**-3    # (M)
gNa = 5

# Concentrations of Potassium (K+)
Co_K = 4*10**-3      # (M)
Ci_K = 150*10**-3    # (M)
gK = 100

# Concentrations of Calcium (Ca2+)
Co_Ca = 1.2*10**-3      # (M)
Ci_Ca = 0.3*10**-6      # (M)

# --- Switch-case style mapping ---
ion_cases = {
    0: ("Sodium",    Co_Na, Ci_Na, z1),
    1: ("Potassium", Co_K,  Ci_K, z1),
    2: ("Calcium",   Co_Ca, Ci_Ca, z2)
}

# --- Function to compute E ---
def compute_E(case):
    """Compute Nernst potential for given ion case."""
    if case not in ion_cases:
        raise ValueError("Invalid case number. Must be 0 or 1.")
    name, Co, Ci, z = ion_cases[case]
    E = ((R * T) / (z * F)) * np.log(Co / Ci)
    return name, E

# --- Compute and store explicitly ---
_, E_Na = compute_E(0)
_, E_K  = compute_E(1)
_, E_Ca = compute_E(2)


# Em at rest
Em = (gNa * E_Na)/(gNa + gK) + (gK * E_K)/(gNa + gK)

# --- Print ion equilibrium potentials ---
print(f"Sodium:       E = {E_Na:.4f} V")
print(f"Potassium:    E = {E_K:.4f} V")
print(f"Calcium:      E = {E_Ca:.4f} V\n")  # <-- newline for spacing
print(f"Resting membraine potential:   E = {Em:.4f} V")
