import math

# --- Given Parameters ---
# Band gap in eV
Eg = 3.0
# 1s exciton resonance peak in eV
E_peak_1s = 1.0
# Target quantum number
n_target = 3

# --- Step 1: Calculate the 1s exciton binding energy ---
Eb_1s = Eg - E_peak_1s

# --- Step 2: Calculate the exciton Rydberg constant (R_ex) for a 2D material ---
# The formula for 2D exciton binding energy is: Eb_n = R_ex / (n - 0.5)^2
# For n=1, we have: Eb_1s = R_ex / (1 - 0.5)^2 = R_ex / 0.25 = 4 * R_ex
# Therefore, R_ex = Eb_1s / 4
R_ex = Eb_1s / 4.0

# --- Step 3: Calculate the binding energy for n=3 ---
# This is the "Rydberg energy for n=3" requested by the user.
Eb_n_target = R_ex / (n_target - 0.5)**2

# --- Final Output ---
# The prompt requires printing the numbers in the final equation.
print(f"The exciton Rydberg constant is R_ex = {R_ex:.4f} eV.")
print(f"The Rydberg energy (binding energy) for n={n_target} is calculated as:")
print(f"E_binding(n={n_target}) = R_ex / ({n_target} - 0.5)^2")
print(f"E_binding(n={n_target}) = {R_ex:.4f} eV / ({n_target} - 0.5)^2 = {Eb_n_target:.4f} eV")