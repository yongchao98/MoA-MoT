import math

# The problem is to count the number of internal adjunctions in a specific 2-category C
# from an object V = F_11^3 to itself.
# As explained in the reasoning, this number is equal to the number of invertible linear maps from V to V,
# which is the size of the general linear group GL(3, F_11).

# Parameters for the calculation
# n is the dimension of the vector space
n = 3
# q is the size of the finite field
q = 11

# The formula for the size of the general linear group GL(n, q) is:
# |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))

# Let's calculate each term for n=3, q=11.
q_n = q**n

# First term: q^n - q^0 = 11^3 - 1
term1 = q_n - q**0
# Second term: q^n - q^1 = 11^3 - 11
term2 = q_n - q**1
# Third term: q^n - q^2 = 11^3 - 121
term3 = q_n - q**2

# The total number of adjunctions is the product of these terms.
result = term1 * term2 * term3

# Print the final equation with the numbers.
print(f"{term1} * {term2} * {term3} = {result}")
