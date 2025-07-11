import sympy as sp

# Define the symbol n for the number of vectors, which goes to infinity.
n = sp.symbols('n')

# As derived in the explanation, the probability P(n) that the sum S satisfies
# ||S||_2 <= sqrt(2) is given by the following expression for large n.
# This expression comes from the Central Limit Theorem applied to the sum of vectors.
# The final equation for the probability is P(n) = c1 - exp(c2/n)
c1 = 1
c2 = -2
Pn = c1 - sp.exp(c2/n)

print(f"The final equation for the probability P(n) for large n is of the form:")
print(f"P(n) = c1 - exp(c2/n)")
print(f"Where the numbers are c1 = {c1} and c2 = {c2}.")
print(f"So, P(n) = {c1} - exp({c2}/n)\n")


# We are asked to find the limit of n * P(n) as n goes to infinity.
expression_to_limit = n * Pn

print("We will compute the limit of the following expression as n -> infinity:")
print(f"Limit(n * ({Pn}))\n")


# Use sympy's limit function to compute the limit of the expression.
limit_value = sp.limit(expression_to_limit, n, sp.oo)

print(f"The value of the limit is:")
# The result from sympy is a symbolic number, which we convert to a standard integer.
print(int(limit_value))