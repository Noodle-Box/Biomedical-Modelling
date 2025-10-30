# --- Physiological constants ---
# Max Saturations (P02)
P02_Rest = 100       # mmHg
P02_Exercise = 125   # mmHg

# Oxygen saturation (P50)
P_healthy = 27
P_sickle = 38


# Case switch system
K_Cases = {
    0: (P02_Rest,      P_healthy),
    1: (P02_Rest,      P_sickle),
    2: (P02_Exercise, P_healthy),
    3: (P02_Exercise, P_sickle),
    }

# --- Switch-Case Function ---
def compute_K(case):
    "Case error"
    if case not in K_Cases:
        raise ValueError("Invalid case number. Must be 0–3.")
    A, B = K_Cases[case]
    p = 2.34
    K = (A**p) / ((B**p) + (A**p))
    return K * 100

K_values = {i: compute_K(i) for i in range(len(K_Cases))}
# --- Compute differences ---
delta_rest      = K_values[0] - K_values[1]  # Rest: healthy - sickle
delta_exercise  = K_values[2] - K_values[3]  # Exercise: healthy - sickle

# --- Printing to terminal ---

for i, K_val in K_values.items():
    print(f"Case {i}: K = {K_val:.1f}%")

print("\n--- Differences ---")
print(f"ΔK (Rest: healthy − sickle)     = {delta_rest:.2f}%")
print(f"ΔK (Exercise: healthy − sickle) = {delta_exercise:.2f}%")