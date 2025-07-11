import math

# Given parameters
Eg = 3.0  # Band gap in eV
E_1s_peak = 1.0  # 1s exciton resonance peak in eV
n = 3     # Target principal quantum number

print("--- Calculation of the Rydberg Energy for the n=3 Exciton State ---")

# Step 1: Calculate the binding energy of the 1s ground state exciton (Eb_1s)
# The relationship is: Eb = Eg - E_peak
Eb_1s = Eg - E_1s_peak
print(f"\nStep 1: Calculate the ground state (n=1) binding energy (Eb_1s).")
print(f"Eb_1s = Band Gap - 1s Resonance Peak")
print(f"Eb_1s = {Eg} eV - {E_1s_peak} eV = {Eb_1s} eV")

# Step 2: Determine the effective 2D Rydberg energy (R_y_2D)
# For a 2D system, the binding energy is Eb_n = R_y_2D / (n - 1/2)^2
# For n=1, this means Eb_1s = R_y_2D / (1 - 0.5)^2, so R_y_2D = Eb_1s * (0.5)^2
# There is a common misconception here. The correct formula relates them as:
# Eb_1s = R_y_2D / (1 - 1/2)^2 = 4 * R_y_2D.
# So, R_y_2D = Eb_1s / 4.
R_y_2D = Eb_1s / (1 / (1 - 0.5)**2)
print(f"\nStep 2: Calculate the effective 2D Rydberg energy (R_y_2D).")
print(f"Using the 2D Rydberg formula for n=1: Eb_1s = R_y_2D / (1 - 0.5)^2")
print(f"Rearranging for R_y_2D gives: R_y_2D = Eb_1s / 4")
print(f"R_y_2D = {Eb_1s} eV / 4 = {R_y_2D} eV")


# Step 3: Calculate the Rydberg energy (binding energy) for the n=3 state (Eb_3)
# We use the 2D Rydberg formula with n=3: Eb_3 = R_y_2D / (3 - 1/2)^2
denominator_n3 = (n - 0.5)**2
Eb_3 = R_y_2D / denominator_n3

print(f"\nStep 3: Calculate the Rydberg energy for n={n}.")
print(f"This is the binding energy of the n={n} state, using the formula Eb_n = R_y_2D / (n - 0.5)^2")
print("\n--- Final Equation ---")
print(f"Eb_{n} = {R_y_2D} eV / ({n} - 0.5)^2")
print(f"Eb_{n} = {R_y_2D} eV / ({n-0.5})^2")
print(f"Eb_{n} = {R_y_2D} eV / {denominator_n3}")
print(f"Eb_{n} = {Eb_3:.3f} eV")
print("------------------------")