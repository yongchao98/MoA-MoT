import math

# Parameters from the problem
# The defect group D is elementary abelian of order 16.
# This means D is a vector space over F_2 of dimension 4.
n = 4  # Dimension of the vector space
q = 2  # Size of the field (from characteristic two)

# The highest possible order of the inertial quotient E is |Aut(D)|,
# which is the order of the General Linear group GL(n, q).

# Calculate the order of GL(n, q) using the formula:
# |GL(n,q)| = (q^n - 1)(q^n - q)(q^n - q^2)...(q^n - q^(n-1))
q_n = q**n
result = 1
terms_symbolic = []
terms_numeric = []

for i in range(n):
    q_i = q**i
    term = q_n - q_i
    result *= term
    terms_symbolic.append(f"({q_n} - {q_i})")
    terms_numeric.append(str(term))

equation_part1 = " * ".join(terms_symbolic)
equation_part2 = " * ".join(terms_numeric)

print("The highest possible order for E is |Aut(D)| = |GL(4, 2)|.")
print("The calculation is as follows:")
print(f"|GL(4, 2)| = {equation_part1}")
print(f"           = {equation_part2}")
print(f"           = {result}")
