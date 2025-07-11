import math

# The problem is to find the value of 1/p_1000.
# From the derivation, the formula for 1/p_n is 4 * cos^2(pi / (n + 2)).

# We set n to 1000.
n = 1000

# The denominator in the angle of the cosine function.
denominator = n + 2

# Calculate the value using the derived formula.
result = 4 * (math.cos(math.pi / denominator))**2

# Print the final equation with all numbers and the calculated result.
# This satisfies the requirement to "output each number in the final equation".
print(f"1/p_{n} = 4 * cos^2(pi / {denominator}) = {result}")