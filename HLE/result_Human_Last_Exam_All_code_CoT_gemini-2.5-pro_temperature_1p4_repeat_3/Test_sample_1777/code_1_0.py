import math

# --- Given Parameters ---
# Band gap in eV
Eg = 3.0
# 1s exciton resonance peak in eV
E_1s = 1.0
# Target quantum number
n_target = 3

# --- Step 1: Calculate the binding energy of the 1s exciton (Eb_1s) ---
Eb_1s = Eg - E_1s
print(f"Step 1: The binding energy of the 1s exciton is Eb_1s = {Eg} eV - {E_1s} eV = {Eb_1s} eV.")

# --- Step 2: Determine the effective Rydberg energy (Ry_eff) for the 2D material ---
# For a 2D system, Eb_n = Ry_eff / (n - 0.5)^2
# For n=1, Eb_1s = Ry_eff / (1 - 0.5)^2 = Ry_eff / 0.25 = 4 * Ry_eff
# Therefore, Ry_eff = Eb_1s / 4
Ry_eff = Eb_1s / 4.0
print(f"Step 2: The effective Rydberg energy is Ry* = {Eb_1s} eV / 4 = {Ry_eff} eV.")

# --- Step 3: Calculate the binding energy for the n=3 state (Eb_3) ---
Eb_3 = Ry_eff / (n_target - 0.5)**2
print(f"\n--- Final Calculation ---")
print(f"The Rydberg energy for n = {n_target} is the binding energy of that state, Eb_{n_target}.")
print(f"The equation is: Eb_n = Ry* / (n - 0.5)^2")
print(f"Plugging in the values: Eb_{n_target} = {Ry_eff:.2f} eV / ({n_target} - 0.5)^2 = {Eb_3:.3f} eV")

<<<0.08>>>