# Define the given physical constants
E_g = 3.0  # Band gap in eV
E_1s = 1.0  # 1s exciton resonance peak in eV
n = 3      # The principal quantum number for the target state

print("This script calculates the exciton binding energy (often called the Rydberg energy) for a specific state in a 2D semiconductor.")
print("-" * 80)

# Step 1: Calculate the binding energy of the 1s exciton (Eb_1s)
Eb_1s = E_g - E_1s
print("Step 1: Calculate the binding energy of the 1s exciton (Eb_1s).")
print(f"   The binding energy is the band gap minus the exciton's resonance energy.")
print(f"   Equation: Eb_1s = E_g - E_1s")
print(f"   Calculation: Eb_1s = {E_g} eV - {E_1s} eV = {Eb_1s} eV\n")

# Step 2: Calculate the effective 2D Rydberg energy (Ry_2D)
# The binding energy for a 2D exciton follows the formula: Eb_n = Ry_2D / (n - 0.5)^2
# For the n=1 state: Eb_1s = Ry_2D / (1 - 0.5)^2 = Ry_2D / (0.5)^2 = 4 * Ry_2D
# We can find Ry_2D by rearranging: Ry_2D = Eb_1s / 4
n_1 = 1
factor_for_n1 = 1 / ((n_1 - 0.5)**2)
Ry_2D = Eb_1s / factor_for_n1
print("Step 2: Calculate the effective 2D Rydberg energy (Ry_2D).")
print(f"   Using the 2D exciton formula for n=1: {Eb_1s} eV = Ry_2D / ({n_1} - 0.5)^2")
print(f"   This simplifies to: {Eb_1s} eV = Ry_2D * {factor_for_n1}")
print(f"   Solving for Ry_2D: Ry_2D = {Eb_1s} eV / {factor_for_n1} = {Ry_2D} eV\n")

# Step 3: Calculate the Rydberg energy for n=3 (i.e., the binding energy Eb_3)
denominator_val = n - 0.5
denominator_sq_val = denominator_val**2
Eb_3 = Ry_2D / denominator_sq_val
print(f"Step 3: Calculate the Rydberg energy for the n={n} state (Eb_{n}).")
print(f"   The formula is: Eb_{n} = Ry_2D / (n - 0.5)^2")
print(f"   Substituting the values into the equation:")
print(f"   Eb_{n} = {Ry_2D} / ({n} - 0.5)^2")
print(f"   Eb_{n} = {Ry_2D} / ({denominator_val})^2")
print(f"   Eb_{n} = {Ry_2D} / {denominator_sq_val}")
print("-" * 80)
print(f"The final calculated Rydberg energy for n={n} is {Eb_3} eV.")
