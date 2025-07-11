import math

n = 1000

# The derived formula for 1/p_n is 4 * cos^2(pi / (n + 2)).
# We need to compute this for n = 1000.
# The numbers in the final equation are 4, pi, and 1002.

four = 4
n_plus_2 = n + 2
pi = math.pi

# Calculate the result
result = four * (math.cos(pi / n_plus_2))**2

# Print the final equation with the numerical values substituted
print(f"1/p_{n} = {four} * cos^2({pi} / {n_plus_2})")

# Print the final calculated value
print(f"The value of 1/p_{n} is: {result}")