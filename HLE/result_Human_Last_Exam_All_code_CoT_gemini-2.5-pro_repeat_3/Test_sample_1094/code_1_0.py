# This script provides the analytical formula for the normalized AC loss
# in an elliptical superconducting bar carrying a transport current.

# The required quantity is the loss per cycle per unit length (Q)
# in the normalized form: 2 * pi * Q / (mu_0 * Ic^2)
# as a function of i = Im / Ic, for i < 1.

# According to the Norris model for an isolated elliptical conductor,
# the formula is:
# Q = (mu_0 * Ic^2 / pi) * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i^2 ]
#
# Rearranging this for the requested normalized form gives:
# 2*pi*Q / (mu_0 * Ic^2) = 2 * [ (1 - i)*ln(1 - i) + (1 + i)*ln(1 + i) - i^2 ]

# The final formula is constructed and printed below.
# The variable 'i' represents the normalized current amplitude Im/Ic.
# The function 'ln' represents the natural logarithm.
# The numbers in the final equation are 2 (overall multiplier), 1 (in two places), and 2 (the exponent).

final_formula = "2 * ((1 - i) * ln(1 - i) + (1 + i) * ln(1 + i) - i**2)"

print("The normalized loss, 2*pi*Q / (mu_0 * Ic^2), as a function of i = Im/Ic is:")
print(final_formula)