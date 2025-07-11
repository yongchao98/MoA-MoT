import math

# This script provides the analytical formula for the normalized AC loss in a
# superconducting elliptic bar for a transport current amplitude below the critical current.

# The problem is to find the normalized loss `2*pi*Q/(mu_0*Ic^2)` as a function of `i = Im/Ic`.
# Based on the Norris model for AC losses in superconductors with sharp edges,
# the final formula is derived as 2 * [(1-i)ln(1-i) + (1+i)ln(1+i) - i^2].

# The following code prints this final expression, showing all the numerical
# constants (2, 1, 1, -1, 2) in the equation as requested.
# The format uses standard mathematical notation where ln is the natural logarithm and i^2 is i squared.

final_formula = "2 * ((1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i**2)"

print("The normalized loss per cycle per unit length, 2*pi*Q/(mu_0*Ic^2), as a function of i = Im/Ic is:")
print(final_formula)
