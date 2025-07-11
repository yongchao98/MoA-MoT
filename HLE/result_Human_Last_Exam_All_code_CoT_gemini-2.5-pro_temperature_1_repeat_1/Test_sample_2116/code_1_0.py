import math

# Based on the analytical derivation, the expected maximum earthquake magnitude
# is given by the expression: pi / (2 * log(2)).
# This script calculates the numerical value of this expression.

# 1. Define the value of pi from the math library.
pi_val = math.pi

# 2. Define the value of the natural logarithm of 2.
log2_val = math.log(2)

# 3. Calculate the denominator of the expression.
denominator = 2 * log2_val

# 4. Calculate the final result.
expected_max_magnitude = pi_val / denominator

# Print the components of the final equation and the result.
print("The final equation for the expected maximum magnitude is: pi / (2 * log(2))")
print("\n--- Equation Components ---")
print(f"Value of pi: {pi_val}")
print(f"Value of log(2): {log2_val}")
print(f"Value of the denominator (2 * log(2)): {denominator}")

print("\n--- Final Result ---")
print(f"Expected Maximum Magnitude = {pi_val} / {denominator}")
print(f"Result: {expected_max_magnitude}")