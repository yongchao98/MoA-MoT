import sys

# Step 1: Define the given physical constants.
E_g = 3.0  # Band gap in eV
E_1s = 1.0  # 1s exciton resonance peak in eV
n1 = 1     # Principal quantum number for the 1s state
n3 = 3     # Principal quantum number for the target state

# Step 2: Calculate the binding energy of the 1s exciton.
# The energy of the resonance peak is E_n = E_g - E_b_n.
# Therefore, the binding energy is E_b_n = E_g - E_n.
E_b_1s = E_g - E_1s

# Step 3: Calculate the effective Rydberg energy (R_y*) for the 2D material.
# The binding energy for a 2D exciton is E_b_n = R_y* / (n - 0.5)^2.
# For n=1: E_b_1s = R_y* / (1 - 0.5)^2
# Therefore, R_y* = E_b_1s * (1 - 0.5)^2
R_y_star = E_b_1s * (n1 - 0.5)**2

# Step 4: Calculate the Rydberg energy (binding energy) for the n=3 state.
# This is the binding energy for the exciton in the n=3 level.
E_b_3s = R_y_star / (n3 - 0.5)**2

# Step 5: Print the results, showing the final equation as requested.
print(f"The binding energy of the 1s exciton is {E_g} eV - {E_1s} eV = {E_b_1s} eV.")
print(f"The effective Rydberg energy R_y* is {E_b_1s} * ({n1} - 0.5)^2 = {R_y_star} eV.")
print("\nThe Rydberg energy for the n=3 state is its binding energy, calculated as:")
print(f"E_b_3s = R_y* / (n - 0.5)^2")
# Final equation with numbers plugged in
print(f"E_b_3s = {R_y_star} / ({n3} - 0.5)^2 = {E_b_3s:.2f} eV")

# Hide the final answer from the output to be placed at the end.
sys.stdout = open('/dev/null', 'w')
print(f'<<<{E_b_3s:.2f}>>>')
