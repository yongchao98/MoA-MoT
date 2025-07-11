import math

# Number of replicas
n_r = 2

# The problem concerns the bosonic sector of the supersymmetric sigma-model
# for symmetry class D with n_r replicas.
# The target manifold is O(2*n_r, 2*n_r) / U(n_r, n_r).
# For n_r = 2, this is O(4,4) / U(2,2).

# The number of non-Grassman variables is the dimension of this manifold.
# dim(G/H) = dim(G) - dim(H)

# Calculate the dimension of G = O(4,4).
# The dimension of O(p,q) is the same as O(p+q).
n_o = 4 + 4
# The dimension of O(n) is n*(n-1)/2.
dim_g = n_o * (n_o - 1) / 2

# Calculate the dimension of H = U(2,2).
# The dimension of U(p,q) is the same as U(p+q).
n_u = 2 + 2
# The dimension of U(n) is n^2.
dim_h = n_u**2

# Calculate the final result
result = dim_g - dim_h

# Print the calculation steps and the final answer.
print(f"The number of non-Grassman variables is the dimension of the manifold O(8)/U(4).")
print(f"Dimension of O(8) = {int(dim_g)}")
print(f"Dimension of U(4) = {int(dim_h)}")
print(f"Total number of variables = {int(dim_g)} - {int(dim_h)} = {int(result)}")
