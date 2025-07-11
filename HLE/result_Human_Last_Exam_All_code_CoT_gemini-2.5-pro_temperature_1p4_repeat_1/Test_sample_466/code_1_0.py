import math

# Set the parameters from the problem
q = 11  # The size of the finite field F_11
n = 3   # The dimension of the vector space F_11^3

# The number of internal adjunctions is equal to the number of
# invertible 3x3 matrices over F_11, which is the order of the
# general linear group GL(n, q).

# The formula for the order of GL(n, q) is:
# |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))

# Calculate the terms of the product
q_n = q**n
term1 = q_n - (q**0)
term2 = q_n - (q**1)
term3 = q_n - (q**2)

# Calculate the final result
result = term1 * term2 * term3

# Print the step-by-step calculation as requested
print("The number of adjunctions is the order of GL(3, F_11).")
print("The calculation is as follows:")
print(f"( {q}^{n} - {q}^0 ) * ( {q}^{n} - {q}^1 ) * ( {q}^{n} - {q}^2 )")
print(f"= ( {q_n} - {q**0} ) * ( {q_n} - {q**1} ) * ( {q_n} - {q**2} )")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")
