import sympy

# Define the variable n for our symbolic calculation
n = sympy.Symbol('n')

# From the analytical derivation, the variance of the sum's components is n/2.
# The probability P(n) is for the norm to be less than or equal to sqrt(2).
# Let R be the radius, so R_squared = 2.
R_squared = 2

# The probability P(n) for large n is approximated by 1 - exp(-R_squared / (2 * variance)).
# The variance is n/2, so the expression becomes 1 - exp(-2 / (2 * n/2)) = 1 - exp(-2/n).
# Let's define the numbers used in the expression for P(n).
constant_in_exp = -2

# The expression for P(n) is P_n = 1 - exp(constant_in_exp / n).
P_n = 1 - sympy.exp(constant_in_exp / n)

# We want to find the limit of n * P(n) as n approaches infinity.
expression_to_limit = n * P_n

# Use sympy to calculate the limit.
limit_value = sympy.limit(expression_to_limit, n, sympy.oo)

# The final equation we are solving is lim_{n->oo} n * (1 - exp(-2/n)).
# The numbers in this equation are 1 (from the subtraction), and -2 in the exponent.
# The result of the limit is the answer.
print(f"The expression for P(n) is approximately: 1 - exp({constant_in_exp}/n)")
print(f"We are calculating the limit of the equation: {expression_to_limit}")
print("The result of the limit is:")
print(limit_value)
