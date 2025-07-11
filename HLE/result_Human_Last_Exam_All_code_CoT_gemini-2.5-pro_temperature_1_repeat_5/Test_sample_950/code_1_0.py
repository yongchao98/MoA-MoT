import math

# This script calculates the rank of the torsion subgroup of the integral
# cohomology ring of the space of 3-subspaces of R^5. This space is
# the Grassmannian Gr(3, 5).

# The argument relies on a theorem from algebraic topology. The Grassmannian
# Gr(k, n) is homotopy equivalent to the homogeneous space G/K where
# G = SO(n) and K = SO(k) x SO(n-k).
# A theorem by A. Borel states that if G and K are compact connected Lie
# groups and rank(G) = rank(K), the integral cohomology ring H*(G/K; Z) is
# torsion-free. A torsion-free ring has a trivial torsion subgroup,
# whose rank is 0.

# We will verify this condition for Gr(3, 5).

# Step 1: Define parameters for the space Gr(k, n).
k = 3
n = 5

def rank_so(m):
    """Calculates the rank of the Lie group SO(m), which is floor(m/2)."""
    return math.floor(m / 2)

# Step 2: Calculate the rank of the group G = SO(n).
rank_G = rank_so(n)
n_minus_k = n - k

print(f"To find the rank of the torsion subgroup of H*(Gr({k}, {n}); Z), we check if the ring is torsion-free.")
print(f"The space Gr({k}, {n}) is homotopy equivalent to the homogeneous space G/K, where G = SO({n}) and K = SO({k}) x SO({n-k}).")
print("According to a theorem by A. Borel, if rank(G) = rank(K), the integral cohomology ring is torsion-free.")
print("-" * 30)

print(f"Step A: Calculate the rank of G = SO({n}).")
# The final equation showing the calculation
print(f"Rank(SO({n})) = floor({n}/2) = {rank_G}")
print("-" * 30)


# Step 3: Calculate the rank of the subgroup K = SO(k) x SO(n-k).
# The rank of a product of groups is the sum of their ranks.
rank_k_group = rank_so(k)
rank_n_minus_k_group = rank_so(n_minus_k)
rank_K = rank_k_group + rank_n_minus_k_group

print(f"Step B: Calculate the rank of K = SO({k}) x SO({n_minus_k}).")
print(f"Rank(SO({k})) = floor({k}/2) = {rank_k_group}")
print(f"Rank(SO({n_minus_k})) = floor({n_minus_k}/2) = {rank_n_minus_k_group}")
# The final equation showing the calculation
print(f"Rank(K) = Rank(SO({k})) + Rank(SO({n_minus_k})) = {rank_k_group} + {rank_n_minus_k_group} = {rank_K}")
print("-" * 30)


# Step 4: Compare the ranks and conclude.
print("Step C: Compare ranks and conclude.")
print(f"Rank(G) = {rank_G}")
print(f"Rank(K) = {rank_K}")

if rank_G == rank_K:
    print("The ranks are equal.")
    print("Therefore, the integral cohomology ring H*(Gr(3, 5); Z) is torsion-free.")
    print("A torsion-free group has a trivial torsion subgroup, {0}, which has a rank of 0.")
    final_answer = 0
else:
    # This path is not taken for the given problem.
    print("The ranks are not equal. This method does not apply.")
    final_answer = "undetermined"

print("\nFinal Answer:")
print(f"The rank of the torsion subgroup is {final_answer}.")
