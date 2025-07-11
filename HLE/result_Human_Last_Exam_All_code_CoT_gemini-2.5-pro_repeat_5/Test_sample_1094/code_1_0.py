# This script presents the formula for the normalized AC loss in an elliptical superconductor.

# The variables in the formula are:
# Q: AC loss per cycle per unit length.
# mu_0: Vacuum permeability.
# Ic: Critical current of the superconductor.
# i: Normalized current amplitude, defined as i = Im/Ic, where Im is the AC current amplitude.
# The condition is i < 1.

# The left-hand side (LHS) of the equation is the standard normalized form for AC loss.
lhs = "2*pi*Q/(mu_0*Ic**2)"

# The right-hand side (RHS) is the resulting function of the normalized current 'i',
# based on the Norris formula for transport current loss in an elliptical wire.
# In mathematical notation, the formula is: 2 * [(1 - i) * ln(1 - i) + i - i^2/2]

# We will construct this formula as a string to be printed.
# The number '2' is the factor resulting from the normalization.
factor = 2

# The terms inside the bracket.
# 'ln' represents the natural logarithm.
# '**' represents the power operator.
term1 = "(1 - i)*ln(1 - i)"
term2 = "i"
term3 = "i**2/2"

# Combine the parts into the final equation string.
final_equation = f"{lhs} = {factor}*({term1} + {term2} - {term3})"

# Print the final equation.
print(final_equation)