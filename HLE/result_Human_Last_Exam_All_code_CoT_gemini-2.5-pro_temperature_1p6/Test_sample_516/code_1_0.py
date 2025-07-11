import math

# The derived formula for the expected value E[T] is 3 - 2 * e^(1/3).
# This script calculates the value of this expression.

# Define the numbers in the final equation
val_a = 3
val_b = -2
val_x = 1/3

# Calculate e^x
exp_val_x = math.exp(val_x)

# Calculate the final result for E[T]
expected_value_T = val_a + val_b * exp_val_x

# Output the final equation and its numerical result
print("The final equation for the expected value of T is derived to be:")
print(f"E[T] = {val_a} + ({val_b}) * e^({val_x})")
print("\nBreaking down the numbers in the equation:")
print(f"The constant term is: {val_a}")
print(f"The coefficient of the exponential term is: {val_b}")
print(f"The exponent is: {val_x}")
print(f"\nThe calculated value is:")
print(f"E[T] = {expected_value_T}")
