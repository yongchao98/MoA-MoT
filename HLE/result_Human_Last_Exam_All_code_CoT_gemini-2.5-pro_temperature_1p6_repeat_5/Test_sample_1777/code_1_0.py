import math

# Define the given physical constants in eV
E_g = 3.0  # Band gap energy
E_1s = 1.0 # 1s exciton resonance peak energy

# Define the quantum numbers for the states
n_1 = 1
n_3 = 3

# Step 1: Calculate the binding energy of the 1s exciton (ground state).
# E_b(1s) = E_g - E_1s
E_b_1s = E_g - E_1s

# Step 2: Calculate the effective Rydberg energy (R_x) for the 2D system.
# In 2D, E_b(n) = R_x / (n - 0.5)^2.
# For n=1, E_b(1s) = R_x / (1 - 0.5)^2 = R_x / 0.25 = 4 * R_x.
# Therefore, R_x = E_b(1s) / 4.
R_x = E_b_1s / ((n_1 - 0.5)**-2)

# Step 3: Calculate the binding energy for the n=3 state.
# This is the "Rydberg energy for n = 3".
# E_b(3) = R_x / (3 - 0.5)^2
denominator_n3 = (n_3 - 0.5)**2
E_b_3 = R_x / denominator_n3

# Print the final result, showing the numbers in the final equation.
print(f"First, the ground state (1s) binding energy is E_b(1s) = {E_g} eV - {E_1s} eV = {E_b_1s} eV.")
print(f"Next, the effective Rydberg energy R_x is found from E_b(1s) = R_x / ({n_1} - 0.5)^2, so R_x = {E_b_1s} eV / 4 = {R_x} eV.")
print("\nFinally, the Rydberg energy for n=3 is calculated:")
print(f"E_b(3) = R_x / (n - 0.5)^2")
print(f"E_b(3) = {R_x:.2f} eV / ({n_3} - 0.5)^2")
print(f"E_b(3) = {R_x:.2f} eV / {denominator_n3}")
print(f"E_b(3) = {E_b_3:.2f} eV")
