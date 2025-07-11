# Based on the analysis of the sequence g_n, the limit g is found
# to have a prime factorization of 2^10 * 3^2 * 5^1.
# This code calculates the final value.

# Exponents
e2 = 10
e3 = 2
e5 = 1

# Base numbers
b2 = 2
b3 = 3
b5 = 5

# Calculate the components
val2 = b2**e2
val3 = b3**e3
val5 = b5**e5

# Calculate the final result
result = val2 * val3 * val5

# Print the final equation with all the numbers
print(f"The limit is {b2}^{e2} * {b3}^{e3} * {b5}^{e5} = {val2} * {val3} * {val5} = {result}")
