import math

# Based on the analytical derivation, the value we need to compute is 3 / (e^2 + 3).
# Let's define the numbers in this final equation.
numerator = 3
exponent_of_e = 2
addend_in_denominator = 3

# Perform the calculation
e_squared = math.exp(exponent_of_e)
denominator_value = e_squared + addend_in_denominator
result = numerator / denominator_value

# The final equation is derived from the expression for -u(0,1)/2.
# We print the final equation including each of its numerical components.
# The equation is: 3 / (e^2 + 3) = result
print(f"The calculation is performed as follows:")
print(f"{numerator} / (e^{exponent_of_e} + {addend_in_denominator})")
print(f"= {numerator} / ({e_squared} + {addend_in_denominator})")
print(f"= {numerator} / {denominator_value}")
print(f"= {result}")