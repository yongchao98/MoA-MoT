import math

# The analytical solution to the problem is E[T] = 3 - 2 * e^(1/3).
# This script calculates the numerical value of this expression.

# Define the constants from the derived formula
constant_term = 3
coefficient = 2
exponent = 1/3

# Calculate the value of e^(1/3)
exp_value = math.exp(exponent)

# Calculate the final expected value
expected_value = constant_term - coefficient * exp_value

# Output the final equation with each number involved
print("The final equation for the expected value E[T] is:")
# The format f-string below prints each number in the final equation
print(f"E[T] = {constant_term} - {coefficient} * {exp_value}")

# Output the final numerical result
print("\nThe computed expected value of T is:")
print(expected_value)