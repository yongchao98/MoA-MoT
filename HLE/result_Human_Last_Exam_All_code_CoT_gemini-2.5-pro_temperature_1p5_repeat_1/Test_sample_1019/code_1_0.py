# The problem asks for the d-threshold for Hamiltonicity, p,
# for a graph H_n with minimum degree d = n/2 - eta.
# The parameters 'n' and 'eta' are symbolic.

# Based on research in probabilistic combinatorics, the threshold p is
# asymptotically equivalent to eta / (2 * n^2).

# The numbers in this final equation are:
# 1: The implicit coefficient of eta in the numerator.
# 2: The coefficient of n^2 in the denominator.
# 2: The exponent of n.

# We will construct and print the formula string, making sure to output
# each of these numbers.

numerator_coefficient = 1
numerator_variable = "eta"
denominator_coefficient = 2
denominator_variable = "n"
denominator_exponent = 2

# The final formula for the d-threshold for Hamiltonicity
final_equation = f"p ~ ({numerator_coefficient} * {numerator_variable}) / ({denominator_coefficient} * {denominator_variable}**{denominator_exponent})"

print("The d-threshold for Hamiltonicity is given by the formula:")
print(final_equation)