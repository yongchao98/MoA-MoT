# Plan:
# The number of internal adjunctions is the order of the general linear group GL(n, F_q)
# for n=3 and q=11. The formula for the order is:
# |GL(n, q)| = (q^n - q^0) * (q^n - q^1) * ... * (q^n - q^(n-1))
# This script will calculate this value.

# n is the dimension of the vector space.
n = 3
# q is the order of the finite field.
q = 11

print(f"The number of internal adjunctions from F_{q}^n to itself is the number of invertible {n}x{n} matrices over F_{q}.")
print(f"This is the order of the general linear group GL({n}, F_{q}).")
print()
print("The calculation for the order of GL(3, 11) is as follows:")

# Calculate the terms of the product formula
q_n = q**n
term1 = q_n - q**0
term2 = q_n - q**1
term3 = q_n - q**2

# Calculate the final result
result = term1 * term2 * term3

# Print the equation with all the numbers, step by step
print(f"|GL({n}, {q})| = ({q}^{n} - {q}^0) * ({q}^{n} - {q}^1) * ({q}^{n} - {q}^2)")
print(f"           = ({q_n} - {q**0}) * ({q_n} - {q**1}) * ({q_n} - {q**2})")
print(f"           = {term1} * {term2} * {term3}")
print(f"           = {result}")