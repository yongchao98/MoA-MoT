import math

# Set the parameters for the field and matrix dimension
q = 11  # The size of the finite field, F_11
n = 3   # The dimension of the vector space, F_11^3

# The number of internal adjunctions is the order of the general linear group GL(n, F_q).
# The formula for the order of GL(n, F_q) is:
# |GL(n, F_q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))

# Calculate each term in the product
q_n = q**n

term1 = q_n - q**0
term2 = q_n - q**1
term3 = q_n - q**2

# Calculate the total number of adjunctions
total_adjunctions = term1 * term2 * term3

# Print the calculation step-by-step
print(f"The number of adjunctions is calculated as the order of GL({n}, F_{q}).")
print(f"The formula is: (q^n - q^0) * (q^n - q^1) * (q^n - q^2)")
print(f"With q={q} and n={n}, q^n = {q}^{n} = {q_n}.")
print(f"The terms are:")
print(f"Term 1: {q_n} - {q**0} = {term1}")
print(f"Term 2: {q_n} - {q**1} = {term2}")
print(f"Term 3: {q_n} - {q**2} = {term3}")
print(f"\nThe final calculation is:")
print(f"{term1} * {term2} * {term3} = {total_adjunctions}")