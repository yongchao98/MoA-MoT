import math

# The problem is solved analytically to find the expression for y(0).
# The final expression for the deflection y at x=0 is y(0) = (3/4) * (9/4)^(3/5).
# This script calculates the numerical value of this expression.

# Define the constants from the derived formula
factor = 3.0 / 4.0
base = 9.0 / 4.0
exponent = 3.0 / 5.0

# Calculate the value of y(0)
y_0 = factor * (base ** exponent)

# Print the final equation with each number and the result
print("The deflection y(0) is calculated using the formula:")
print(f"y(0) = ({factor}) * ({base})^({exponent})")
print("\nThe calculated deflection at x = 0 is:")
print(f"y(0) = {y_0}")