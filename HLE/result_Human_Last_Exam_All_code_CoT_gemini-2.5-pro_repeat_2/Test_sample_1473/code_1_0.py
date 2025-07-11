import math

# The value of the integral I is given by the expression pi * ln(1 + sqrt(2))
# This script calculates and prints this value, showing the components of the final equation.

# Define the constants and intermediate values
pi_val = math.pi
sqrt2_val = math.sqrt(2)
term1 = 1 + sqrt2_val
term2 = math.log(term1) # math.log is the natural logarithm (ln)

# Calculate the final result
I = pi_val * term2

# Print the breakdown of the final equation
print("The integral I is calculated as: I = pi * ln(1 + sqrt(2))")
print("\n--- Values of each component ---")
print(f"pi = {pi_val}")
print(f"sqrt(2) = {sqrt2_val}")
print(f"1 + sqrt(2) = {term1}")
print(f"ln(1 + sqrt(2)) = {term2}")
print("\n--- Final Result ---")
print(f"I = {pi_val} * {term2}")
print(f"I = {I}")
