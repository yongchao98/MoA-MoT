import math

# The problem asks to find 1/p_n for n=1000.
# From the derivation, we found the general formula for p_n is:
# p_n = 1 / (4 * cos^2(pi / (n + 2)))
# Therefore, 1/p_n = 4 * cos^2(pi / (n + 2))

# Set the value of n
n = 1000

# Calculate the denominator for the angle
denominator = n + 2

# Calculate the angle in radians
angle = math.pi / denominator

# Calculate the value of 1/p_1000
result = 4 * (math.cos(angle))**2

# As requested, output the numbers in the final equation leading to the result.
# The final equation is 1/p_1000 = 4 * cos^2(pi / 1002)
print(f"The equation for 1/p_{n} is: 4 * cos^2(pi / {denominator})")
print(f"The value of 1/p_{n} for n = {n} is: {result}")