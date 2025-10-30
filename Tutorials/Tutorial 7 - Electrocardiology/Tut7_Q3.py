import numpy as np
import matplotlib.pyplot as plt

## Parameters 

F = 9.649*10**4     # (C/mol)
R = 8.314           # (J/mol*k)
z = 1
T = 310             # (Kelvin, Temp)

# Conductance graphs:
Data = np.array([[0.01, 0.4, 0.8, 1.2, 1.6, 2], 
                 [0.366751, 2.186175, 5.478086, 9.486194, 13.48798, 17.05401], 
                 [0.010589, 37.48381, 31.23947, 21.33964, 14.36902, 9.674778]])


# Concentrations of Sodium
Co_Na = 145*10**-3   # (M)
Ci_Na = 12*10**-3    # (M)

# Concentrations of Potassium
Co_K = 4*10**-3      # (M)
Ci_K = 155*10**-3    # (M)


# --- Switch-case style mapping ---
ion_cases = {
    0: ("Sodium",    Co_Na, Ci_Na),
    1: ("Potassium", Co_K,  Ci_K),
}

# --- Function to compute E ---
def compute_E(case):
    """Compute Nernst potential for given ion case."""
    if case not in ion_cases:
        raise ValueError("Invalid case number. Must be 0 or 1.")
    name, Co, Ci = ion_cases[case]
    E = ((R * T) / (z * F)) * np.log(Co / Ci)
    return name, E

# --- Compute and store explicitly ---
_, E_Na = compute_E(0)
_, E_K  = compute_E(1)


# --- Separate rows ---
time = Data[0]
gK = Data[1]
gNa  = Data[2]


# --- Compute Em iteratively and save to a list ---
Em_list = []

for t in range(len(time)):
    
    E_t = (gNa[t] * E_Na)/(gNa[t] + gK[t]) + (gK[t] * E_K)/(gNa[t] + gK[t])
    Em_list.append(E_t)

# --- Print ion equilibrium potentials ---
print(f"Sodium:    E = {E_Na:.3f} V")
print(f"Potassium: E = {E_K:.3f} V\n")  # <-- newline for spacing

# --- Header for time-series results ---
print("Time (s)   Membrane Potential (V)")
print("-" * 35)

# --- Print membrane potential over time ---
for t, E in zip(time, Em_list):
    print(f"{t:>6.1f}      {E:>8.4f}")


## Plotting ##
plt.figure()
plt.subplot(211)
plt.plot(time, gK, marker='o', label='gK (K⁺)', color='b')
plt.plot(time, gNa, marker='s', label='gNa (Na⁺)', color='r')
plt.ylabel('Conductance (mS/cm^2)')
plt.subplot(212)
plt.plot(time, Em_list, marker='d', color='k')
plt.xlabel('Time (s)')
plt.ylabel('Membrane Potential (V)')
plt.show()