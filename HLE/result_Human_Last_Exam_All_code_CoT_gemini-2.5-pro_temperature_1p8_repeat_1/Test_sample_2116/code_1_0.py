import math

# The problem boils down to calculating the expression pi / (2 * log(2)).
# Here, we will compute this value and show the components of the final equation.

# Get the value of pi and natural logarithm of 2 from the math library.
pi_val = math.pi
log2_val = math.log(2)

# Calculate the numerator and denominator of the final expression.
numerator = pi_val
denominator = 2 * log2_val

# Calculate the expected maximum magnitude.
expected_max_magnitude = numerator / denominator

# Print the final equation with the computed values.
print(f"The expected maximum earthquake magnitude is given by the equation:")
print(f"E[Max Magnitude] = \u03C0 / (2 * log(2))")
print(f"E[Max Magnitude] = {numerator} / (2 * {log2_val})")
print(f"E[Max Magnitude] = {numerator} / {denominator}")
print(f"E[Max Magnitude] = {expected_max_magnitude}")
