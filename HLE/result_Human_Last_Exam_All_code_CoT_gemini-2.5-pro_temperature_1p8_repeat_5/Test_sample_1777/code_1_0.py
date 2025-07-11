import math

# --- Given Parameters ---
# Band gap in eV
E_g = 3.0
# Resonance peak for the 1s exciton in eV
E_1s = 1.0
# Principal quantum number for the target state
n_target = 3

# --- Step 1: Calculate the binding energy of the 1s state ---
# The binding energy is the difference between the band gap and the exciton energy level.
E_b_1s = E_g - E_1s
print(f"The binding energy of the 1s state (E_b_1s) is calculated as:")
print(f"E_b_1s = E_g - E_1s = {E_g} eV - {E_1s} eV = {E_b_1s} eV\n")

# --- Step 2: Calculate the 2D Rydberg energy (R_X) ---
# In 2D materials, the exciton binding energy follows the formula: E_b_n = R_X / (n - 0.5)^2
# We can find R_X using the 1s state (n=1).
n_1 = 1
denominator_1s = (n_1 - 0.5)**2
R_X = E_b_1s * denominator_1s
print(f"The 2D Rydberg Energy (R_X) is found using the 1s state:")
print(f"R_X = E_b_1s * (n - 0.5)^2 = {E_b_1s} eV * ({n_1} - 0.5)^2 = {R_X} eV\n")

# --- Step 3: Calculate the Rydberg energy (binding energy) for n = 3 ---
# Now, we use the calculated R_X to find the binding energy for the n=3 state.
denominator_3s = (n_target - 0.5)**2
E_b_3s = R_X / denominator_3s
print(f"The Rydberg energy for n = {n_target} is calculated as:")
print(f"E_b_{n_target} = R_X / (n - 0.5)^2")
print(f"E_b_{n_target} = {R_X} eV / ({n_target} - 0.5)^2")
print(f"E_b_{n_target} = {R_X} eV / {denominator_3s}")
print(f"Final Answer: The Rydberg energy for n = {n_target} is {E_b_3s:.3f} eV.")
