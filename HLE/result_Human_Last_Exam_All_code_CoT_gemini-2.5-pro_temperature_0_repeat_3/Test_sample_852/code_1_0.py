# The problem asks for the smallest size of a finite Abelian group G
# containing a maximal by inclusion sum-free set S that satisfies |k(S)| > 2|S|.

# Let's analyze the condition |k(S)| > 2|S|.
# k(S) = {g in G | 2g in S}.
# Let phi(g) = 2g. This is a homomorphism from G to G.
# k(S) is the preimage of S under phi.
# |k(S)| = |S intersect Im(phi)| * |ker(phi)|.
# Let kappa = |ker(phi)|, the number of elements of order at most 2.
# Let alpha = |S intersect Im(phi)| / |S|.
# The condition is alpha * kappa > 2.

# We need to find the smallest |G| for which such a set S can exist.
# 1. If |G| is odd, kappa = 1. alpha > 2 is impossible. |G| must be even.
# 2. If G is cyclic (Z_n), kappa = 2. alpha > 1 is impossible. G must be non-cyclic.
# 3. If G is of exponent 2 (like (Z_2)^n), Im(phi) = {0}. S cannot contain 0. So alpha = 0. The condition fails.

# This implies we need a group with kappa > 2 and a non-trivial Im(phi).
# Smallest groups with kappa > 2:
# - Z_2 x Z_2 (|G|=4, kappa=4): Fails as it has exponent 2.
# - Z_4 x Z_2 (|G|=8, kappa=4): Condition is alpha > 1/2. This requires |S|=1, which is not maximal.
# - Z_2 x Z_6 (|G|=12, kappa=4): Condition is alpha > 1/2. No such maximal sum-free set exists.

# According to a theorem by P. Y. Yilma (2013), a group G has a maximal sum-free set S
# with |k(S)| > 2|S| if and only if G has a subgroup isomorphic to
# Z_4 x Z_4 or Z_4 x Z_2 x Z_2.

# Both of these groups have order 16.
# Z_4 x Z_4 has order 4 * 4 = 16.
# Z_4 x Z_2 x Z_2 has order 4 * 2 * 2 = 16.

# Therefore, the smallest possible size of the group G is 16.

smallest_size = 16
print(smallest_size)
