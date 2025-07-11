import math

q = 11
n = 3

# Calculate the terms in the formula for the order of GL(n, q)
q_n = q**n
term1 = q_n - q**0
term2 = q_n - q**1
term3 = q_n - q**2

# Calculate the total number of adjunctions
result = term1 * term2 * term3

# Print the equation with all numbers
print(f"The number of internal adjunctions is the order of GL({n}, F_{q}), calculated as:")
print(f"({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2) = ({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
print(f"= {term1} * {term2} * {term3}")
print(f"= {result}")
