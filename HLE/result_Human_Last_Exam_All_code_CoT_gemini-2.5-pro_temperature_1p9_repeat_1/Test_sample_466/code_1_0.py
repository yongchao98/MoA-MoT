# Parameters for the calculation
q = 11
n = 3

# Calculate the size of GL(n, q)
# The formula is the product of (q^n - q^i) for i from 0 to n-1.

# First term: (q^n - q^0) = (11^3 - 1)
q_n = q**n
term1 = q_n - (q**0)

# Second term: (q^n - q^1) = (11^3 - 11)
term2 = q_n - (q**1)

# Third term: (q^n - q^2) = (11^3 - 121)
term3 = q_n - (q**2)

# Calculate the final result
result = term1 * term2 * term3

# Print the calculation step by step
print(f"The number of adjunctions is the size of GL({n}, F_{q}), which is calculated as:")
print(f"(q^n - q^0) * (q^n - q^1) * (q^n - q^2)")
print(f"= ({q}^{n} - {q**0}) * ({q}^{n} - {q**1}) * ({q}^{n} - {q**2})")
print(f"= ({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")
