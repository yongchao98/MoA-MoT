import math

# The problem asks for the constant 'b' in the asymptotic formula C(n) ~ b * n^(3/2)
# for the expected time of a certain random walk process on a random n-vertex tree.

# Based on analyzing the scaling laws of various quantities related to random walks on random trees,
# the quantity C(n) that scales as n^(3/2) is the expected commute time between two randomly
# chosen vertices, averaged over all random trees.

# The derivation is as follows:
# 1. The commute time K(u,v) on a tree with n vertices is 2*(n-1)*d(u,v), where d(u,v) is the distance.
# 2. The expected distance E[d(u,v)] between two random vertices on a random tree is ~ sqrt(pi*n/2).
# 3. Thus, the expected commute time E[K(u,v)] ~ 2*n * E[d(u,v)] ~ 2*n * sqrt(pi*n/2) = sqrt(4*n^2 * pi*n/2) = sqrt(2*pi*n^3) = sqrt(2*pi) * n^(3/2).
# 4. Therefore, the constant b is sqrt(2*pi).

b_squared = 2 * math.pi
b = math.sqrt(b_squared)

print("The asymptotic behavior is C(n) ~ b * n^(3/2)")
print(f"The exact value of the constant b is sqrt(2 * pi)")
print(f"Let pi = {math.pi}")
print(f"b = sqrt(2 * {math.pi})")
print(f"b = sqrt({b_squared})")
print(f"b = {b}")
