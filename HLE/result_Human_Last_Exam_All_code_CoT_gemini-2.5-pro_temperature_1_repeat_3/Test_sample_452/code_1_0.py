import math

# The problem is to find the constant b in the asymptotic formula C(n) ~ b * n^(3/2)
# for the expected cover and return time of a random walk on a uniform random tree.

# Based on the derivation, the constant b is equal to sqrt(pi / 8).
# Here we calculate this value.

pi = math.pi
b_squared = pi / 8
b = math.sqrt(b_squared)

print("The constant b is the square root of (pi / 8).")
print(f"pi = {pi}")
print(f"b = sqrt({pi} / 8)")
print(f"b = {b}")
