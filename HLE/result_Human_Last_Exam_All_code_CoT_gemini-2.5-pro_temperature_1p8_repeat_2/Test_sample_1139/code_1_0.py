import math

# Number of replicas
N = 2

# The number of non-Grassman (bosonic) variables for the supersymmetric sigma-model
# of symmetry class D is given by the dimension of the corresponding bosonic symmetric space.
# For N replicas, this space is O(2*N) / (O(N) x O(N)).
# The dimension is calculated as: dim(O(2*N)) - dim(O(N) x O(N)).

# The dimension of the orthogonal group O(n) is n*(n-1)/2.

# Dimension of the group G = O(2*N)
group_G_dim = (2 * N) * (2 * N - 1) // 2

# Dimension of the subgroup H = O(N) x O(N)
# The dimension of a direct product of groups is the sum of their dimensions.
# dim(H) = dim(O(N)) + dim(O(N)) = 2 * dim(O(N))
subgroup_H_dim = 2 * (N * (N - 1) // 2)

# The result is the dimension of the quotient space G/H.
result = group_G_dim - subgroup_H_dim

# --- Output the result as requested ---
print(f"The problem is to find the number of non-Grassman variables for a model with N={N} replicas in symmetry class D.")
print("This number is the dimension of the symmetric space O(2N)/[O(N)xO(N)].")
print(f"The calculation is: dim(O({2*N})) - (dim(O({N})) + dim(O({N})))")
print("\nFirst, we calculate the dimension of the individual groups:")
print(f"dim(O({2*N})) = ({2*N}*({2*N}-1))/2 = {group_G_dim}")
print(f"dim(O({N})) = ({N}*({N}-1))/2 = {N * (N - 1) // 2}")
print(f"So, dim(O({N})xO({N})) = {subgroup_H_dim}")
print("\nFinally, we subtract the subgroup dimension from the group dimension.")
print("The final equation is:")
print(f"{group_G_dim} - {subgroup_H_dim} = {result}")
