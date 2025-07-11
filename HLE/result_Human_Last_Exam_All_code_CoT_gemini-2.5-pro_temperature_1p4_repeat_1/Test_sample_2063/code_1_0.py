# Plan: The limiting probability is given by (7/8)^8.
# This code calculates this value and prints the components of the expression.

# The base of the fraction
base = 7
# The divisor of the fraction
divisor = 8
# The power to which the fraction is raised
power = 8

# Calculate the numerator and denominator
numerator = base**power
denominator = divisor**power

# Calculate the final result
result = numerator / denominator

# Print the final equation with all its numbers
print(f"The limiting probability is P(n) as n -> infinity.")
print(f"P = ({base}/{divisor})^{power}")
print(f"P = {numerator} / {denominator}")
print(f"P = {result}")
