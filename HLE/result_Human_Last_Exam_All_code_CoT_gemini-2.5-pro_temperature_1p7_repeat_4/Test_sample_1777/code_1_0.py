import math

# --- Given values ---
E_g = 3.0  # Band gap in eV
E_1s = 1.0  # 1s exciton resonance peak in eV
n = 3      # Principal quantum number for the target state

# Step 1: Calculate the ground state (n=1) binding energy
E_b_1s = E_g - E_1s
print(f"Step 1: Calculate the ground state (n=1) binding energy.")
print(f"E_b(n=1) = E_g - E_1s = {E_g} eV - {E_1s} eV = {E_b_1s} eV")
print("-" * 30)

# Step 2: Determine the effective 2D Rydberg energy (R_y_2D)
# For a 2D exciton, E_b(n=1) = R_y_2D / (1 - 0.5)^2 = 4 * R_y_2D
R_y_2D = E_b_1s / 4.0
print(f"Step 2: Determine the effective 2D Rydberg energy (R_y_2D).")
print(f"R_y_2D = E_b(n=1) / 4 = {E_b_1s} eV / 4 = {R_y_2D} eV")
print("-" * 30)

# Step 3: Calculate the binding energy for n=3
# This is the "Rydberg energy for n=3"
denominator = n - 0.5
E_b_n3 = R_y_2D / (denominator**2)
print(f"Step 3: Calculate the binding energy for n={n}.")
print("The formula is: E_b(n) = R_y_2D / (n - 0.5)^2")
print(f"E_b(n={n}) = {R_y_2D} eV / ({n} - 0.5)^2")
print(f"        = {R_y_2D} eV / ({denominator})^2")
print(f"        = {R_y_2D} eV / {denominator**2}")
print(f"        = {E_b_n3} eV")
print("-" * 30)

print(f"The final Rydberg energy for n = {n} is: {E_b_n3} eV")