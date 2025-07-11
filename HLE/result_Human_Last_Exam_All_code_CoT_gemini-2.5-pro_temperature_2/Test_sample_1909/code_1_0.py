import math

# The problem asks for the limit of the expected ratio of remaining numbers to n,
# as n approaches infinity.
# Let E_n be the expected number of remaining numbers.
# The limit is given by the constant A in the asymptotic expression E_n ~ A*(n+2).
# This constant A is found by evaluating the integral of (1-x)^2 * exp(2x) from 0 to 1.
# The value of this integral is (e^2 - 5) / 4.

# Let's calculate this value.

# Define the constants from the final equation
e_val = math.e
power = 2
subtrahend = 5
divisor = 4

# Perform the calculation
e_squared = e_val ** power
numerator = e_squared - subtrahend
result = numerator / divisor

# Print the components and the final result
print("The limit is given by the expression: (e^2 - 5) / 4")
print(f"Value of e: {e_val}")
print(f"Value of the power: {power}")
print(f"Value of the subtrahend: {subtrahend}")
print(f"Value of the divisor: {divisor}")
print(f"Intermediate calculation e^2: {e_squared}")
print(f"Intermediate calculation e^2 - 5: {numerator}")
print(f"Final result (e^2 - 5) / 4: {result}")