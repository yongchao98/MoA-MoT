import math

# --- Given parameters ---
# Band gap of the 2D semiconductor in eV
E_g = 3.0
# Resonance peak energy for the 1s exciton in eV
E_1s = 1.0
# The principal quantum number for the target state
n_target = 3

# --- Step 1: Calculate the binding energy of the 1s exciton ---
# The formula is: E_binding = E_g - E_exciton
E_b_1s = E_g - E_1s

print("Step 1: Calculate the 1s exciton binding energy (E_b,1s)")
print(f"E_b,1s = E_g - E_1s")
print(f"E_b,1s = {E_g} eV - {E_1s} eV = {E_b_1s} eV")
print("-" * 30)

# --- Step 2: Calculate the effective Rydberg energy (R_y) for the exciton ---
# For a 2D exciton, the binding energy formula is: E_b,n = R_y / (n - 1/2)^2
# For n=1, this simplifies to E_b,1s = R_y / (1 - 0.5)^2 = R_y / 0.25 = 4 * R_y
# Therefore, R_y = E_b,1s / 4
R_y = E_b_1s / 4.0

print("Step 2: Calculate the effective Rydberg energy (R_y)")
print("From the 2D exciton model, E_b,1s = R_y / (1 - 0.5)^2 = 4 * R_y")
print(f"R_y = E_b,1s / 4")
print(f"R_y = {E_b_1s} eV / 4 = {R_y} eV")
print("-" * 30)

# --- Step 3: Calculate the Rydberg energy (binding energy) for n=3 ---
# Use the 2D exciton binding energy formula with n = 3 and the calculated R_y
denominator = (n_target - 0.5)**2
E_b_n_target = R_y / denominator

print(f"Step 3: Calculate the Rydberg energy (binding energy) for n={n_target} (E_b,{n_target})")
print(f"E_b,n = R_y / (n - 1/2)^2")
print(f"E_b,{n_target} = {R_y} eV / ({n_target} - 0.5)^2")
print(f"E_b,{n_target} = {R_y} eV / ({n_target - 0.5})^2")
print(f"E_b,{n_target} = {R_y} eV / {denominator}")
print(f"E_b,{n_target} = {E_b_n_target:.3f} eV")
print("-" * 30)

print(f"The final answer is {E_b_n_target:.3f} eV.")
