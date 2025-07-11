import sympy

# Define the variable for the total number of vectors
n = sympy.Symbol('n')

# From the analysis based on the Central Limit Theorem, we derived the
# approximate probability P(n) for large n.
# The variance of the components of the sum vector S is sigma_sq = n / 2.
# The probability P(n) = P(||S||^2 <= 2) is given by P(n) ~ 1 - exp(-1 / sigma_sq).
# Substituting sigma_sq = n/2 gives P(n) ~ 1 - exp(-2/n).

# We want to find the limit of the expression n * P(n).
# The constant 'c' in the exponent comes from the norm threshold and the variance calculation.
# c = (threshold_r^2) / (variance_coefficient) = 2 / (1/2) = 4 is wrong.
# Let's check my derivation again. P(n) = 1-exp(-1/sigma^2) = 1 - exp(-1/(n/2)) = 1-exp(-2/n).
# So the constant in the numerator is 2.
c = 2

# Define the expression for P(n)
P_n = 1 - sympy.exp(-c / n)

# Define the full expression whose limit is required
expression_to_limit = n * P_n

# Now, we compute the limit as n -> infinity.
limit_value = sympy.limit(expression_to_limit, n, sympy.oo)

print("The problem asks for the limit of n * P(n) as n tends to infinity.")
print("Using a Central Limit Theorem approximation, the expression is:")
# The print function will display the symbolic expression
print(f"lim_{n->oo} {expression_to_limit}")
# The following line shows the numbers in the final equation.
print(f"Here, the constant in the exponent of the equation for P(n) is {c}.")
print("\nThe calculated value of the limit is:")
print(limit_value)
