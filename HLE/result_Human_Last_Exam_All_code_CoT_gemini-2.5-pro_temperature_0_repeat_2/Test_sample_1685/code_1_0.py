import math

# Based on the analysis, the solution for large x is approximated by y(x) = A/x,
# where the coefficient A is determined by the equation A^3 = -6.

# The constant in the equation is -6.
constant_value = -6.0

# We need to calculate A = (-6)^(1/3).
coefficient_A = -math.pow(abs(constant_value), 1/3)

# Round the coefficient to two decimal places as requested.
rounded_coefficient_A = round(coefficient_A, 2)

# The final analytical expression is y(x) = A / x.
# We will print the equation with the calculated and rounded coefficient.
print("The approximate analytical expression for y(x) in the large x regime is:")
print(f"y(x) = {rounded_coefficient_A:.2f} / x")
