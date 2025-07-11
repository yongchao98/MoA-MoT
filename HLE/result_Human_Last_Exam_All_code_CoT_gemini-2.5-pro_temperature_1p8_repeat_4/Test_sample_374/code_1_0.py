# Parameters from the problem
# The defect group D is elementary abelian of order 16.
# This corresponds to a vector space of dimension n=4 over the field F_q with q=2 elements.
n = 4
q = 2

# The highest possible order for the inertial quotient E is the order of Aut(D),
# which is isomorphic to the order of the general linear group GL(n, q).

# Calculate the terms of the formula for the order of GL(n, q):
# |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))
term1 = q**n - q**0
term2 = q**n - q**1
term3 = q**n - q**2
term4 = q**n - q**3

# Calculate the final product
highest_order = term1 * term2 * term3 * term4

# Print the final equation with all its components as requested.
print("The highest possible order for E is |GL(4, 2)|.")
print("The calculation is:")
print(f"{term1} * {term2} * {term3} * {term4} = {highest_order}")