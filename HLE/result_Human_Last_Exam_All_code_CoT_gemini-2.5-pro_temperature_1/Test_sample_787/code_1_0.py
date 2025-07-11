# Based on the mathematical derivation, the limit of the sequence is 2^10 * 3^2 * 5^1.
# This script calculates this value and prints it in the required format.

# Exponents determined from the analysis
A = 10
B = 2
C = 1

# Calculate the components
power_of_2 = 2**A
power_of_3 = 3**B
power_of_5 = 5**C

# Calculate the final limit
limit_g_n = power_of_2 * power_of_3 * power_of_5

# Print the result in the specified equation format
print(f"{power_of_2} * {power_of_3} * {power_of_5} = {limit_g_n}")
print("The equation representing the limit is:")
print(f"2**{A} * 3**{B} * 5**{C} = {limit_g_n}")
