import sys

# Define the given physical parameters
Eg = 3.0  # Band gap in eV
E1 = 1.0  # 1s exciton resonance peak energy in eV
n = 3     # Target principal quantum number

# Step 1: Calculate the binding energy of the 1s exciton (Eb_1)
# The binding energy is the difference between the band gap and the exciton resonance energy.
# Eb_1 = Eg - E1
Eb_1 = Eg - E1

print(f"Step 1: Calculate the binding energy of the 1s exciton (Eb_1).")
print(f"The formula is: Eb_1 = Eg - E1")
print(f"Eb_1 = {Eg} eV - {E1} eV = {Eb_1} eV\n")

# Step 2: Determine the effective 2D Rydberg energy (Ry_2D)
# For a 2D exciton, the binding energy series is Eb_n = Ry_2D / (n - 0.5)^2.
# For the ground state (n=1), this gives Eb_1 = Ry_2D / (1 - 0.5)^2 = 4 * Ry_2D.
# We can find Ry_2D from Eb_1.
Ry_2D = Eb_1 / 4.0

print(f"Step 2: Determine the effective 2D Rydberg energy (Ry_2D).")
print(f"From the 2D exciton model, Ry_2D = Eb_1 / (1 - 0.5)^-2 = Eb_1 / 4.")
print(f"Ry_2D = {Eb_1} eV / 4 = {Ry_2D} eV\n")

# Step 3: Calculate the Rydberg energy for n=3 (Eb_3)
# This is the binding energy of the n=3 exciton state.
denominator = (n - 0.5)**2
Eb_3 = Ry_2D / denominator

print(f"Step 3: Calculate the Rydberg energy for n = {n} (Eb_{n}).")
print(f"The formula is: Eb_n = Ry_2D / (n - 0.5)^2")
# Using f-strings to format the output string with the variable values.
# Each number in the final equation is printed.
equation = f"Eb_{n} = {Ry_2D} eV / ({n} - 0.5)^2 = {Ry_2D} eV / {denominator} = {Eb_3:.2f} eV"
print(equation)

# Final answer in the required format
# We use sys.stdout.write to avoid adding an extra newline
sys.stdout.write(f"\n<<<{Eb_3:.2f}>>>")