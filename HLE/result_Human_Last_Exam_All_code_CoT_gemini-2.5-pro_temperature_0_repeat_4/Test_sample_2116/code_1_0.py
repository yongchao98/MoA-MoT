import math

# The problem of finding the expected maximum earthquake magnitude under these
# conditions simplifies to the elegant mathematical expression: pi / (2 * log(2)).
# Here, we will calculate this value.

# Define the constants from the equation
pi = math.pi
two = 2
log_of_2 = math.log(2)

# Calculate the final result
expected_max_magnitude = pi / (two * log_of_2)

# As requested, here is the final equation with each number printed out,
# followed by the final result.
print("Final Equation: pi / (2 * log(2))")
print(f"Calculation: {pi} / ({two} * {log_of_2})")
print(f"Result: {expected_max_magnitude}")