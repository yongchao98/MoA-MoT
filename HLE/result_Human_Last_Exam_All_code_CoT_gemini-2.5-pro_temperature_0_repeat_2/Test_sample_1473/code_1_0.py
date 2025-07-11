import math

# The integral I has been analytically solved, and the result is pi * ln(1 + sqrt(2)).
# This script calculates the numerical value of this expression.

# The final equation for the integral's value is I = pi * ln(1 + sqrt(2)).
# We will print the values of the components of this equation.

pi_val = math.pi
one_val = 1.0
sqrt2_val = math.sqrt(2)

# Argument of the natural logarithm
ln_argument = one_val + sqrt2_val

# Value of the natural logarithm term
ln_val = math.log(ln_argument)

# The final value of the integral I
integral_value = pi_val * ln_val

print("The final equation is I = pi * ln(1 + sqrt(2))")
print(f"Value of pi: {pi_val}")
print(f"Value of 1: {one_val}")
print(f"Value of sqrt(2): {sqrt2_val}")
print(f"Value of the term (1 + sqrt(2)): {ln_argument}")
print(f"Value of ln(1 + sqrt(2)): {ln_val}")
print(f"The final value of the integral I is: {integral_value}")
