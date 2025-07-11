import math

# Set the parameters for the finite field and vector space dimension
q = 11  # The size of the finite field F_q
n = 3   # The dimension of the vector space F_q^n

# Calculate the terms in the product formula for the size of GL(n, q)
q_n = q**n
term1 = q_n - q**0
term2 = q_n - q**1
term3 = q_n - q**2

# Calculate the final result
result = term1 * term2 * term3

# Print the explanation and the step-by-step calculation
print("An adjunction is uniquely determined by an invertible linear map f from F_11^3 to itself.")
print("The number of such maps is the size of the general linear group GL(3, 11).")
print(f"The formula for the size of GL(n, q) is (q^n - q^0) * ... * (q^n - q^(n-1)).")
print(f"For n={n} and q={q}, the calculation is:")
print(f"|GL({n}, {q})| = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2)")
print(f"           = ({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
print(f"           = {term1} * {term2} * {term3}")
print(f"           = {result}")
