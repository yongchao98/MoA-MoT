import math

# --- Given parameters ---
# Band gap in eV
Eg = 3.0
# 1s exciton resonance peak energy in eV
E_1s = 1.0
# Target principal quantum number
n_target = 3

# --- Step 1: Calculate the exciton Rydberg energy (R_y) ---
# The energy of the 1s exciton (n=1) is given by:
# E_1s = Eg - R_y / (1 - 0.5)^2  =>  E_1s = Eg - 4 * R_y
# We can rearrange to solve for R_y:
# 4 * R_y = Eg - E_1s
# R_y = (Eg - E_1s) / 4
R_y = (Eg - E_1s) / 4.0

# --- Step 2: Calculate the binding energy for n=3 ---
# The "Rydberg energy for n=3" is the binding energy of the n=3 state.
# E_binding(n) = R_y / (n - 0.5)^2
denominator = (n_target - 0.5)**2
E_binding_3 = R_y / denominator

# --- Print the results ---
print("The problem is to find the binding energy of the n=3 exciton state.")
print("\nStep 1: Calculate the exciton Rydberg energy, R_y.")
print(f"From the 1s exciton energy, we have the equation: {E_1s} eV = {Eg} eV - 4 * R_y")
print(f"Solving for R_y: R_y = ({Eg} eV - {E_1s} eV) / 4 = {R_y} eV")

print("\nStep 2: Calculate the binding energy for n = 3, E_binding(3).")
print(f"The formula is: E_binding(n) = R_y / (n - 0.5)^2")
print(f"Plugging in the values for n = {n_target}:")
print(f"E_binding({n_target}) = {R_y} eV / ({n_target} - 0.5)^2")
print(f"E_binding({n_target}) = {R_y} eV / {denominator}")
print(f"The Rydberg energy for n = 3 is: {E_binding_3:.3f} eV")