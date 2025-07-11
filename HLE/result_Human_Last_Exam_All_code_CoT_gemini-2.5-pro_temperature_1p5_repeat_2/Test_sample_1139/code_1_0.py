import math

# Number of replicas
N_r = 2
print(f"The number of replicas is N_r = {N_r}")

# The bosonic sector of the supersymmetric sigma-model for class D with N_r replicas
# is parametrized by the manifold M = [SO(4*N_r)/U(2*N_r)] x [SO(4*N_r)/U(2*N_r)].
# We need to calculate its dimension.

# Step 1: Determine the parameters of the Lie groups.
m = 4 * N_r
k = 2 * N_r
print(f"The manifold for one component is SO(m)/U(k) with m = 4*N_r = {m} and k = 2*N_r = {k}.")
print("-" * 20)

# Step 2: Calculate the dimension of the special orthogonal group SO(m).
# dim(SO(m)) = m * (m - 1) / 2
dim_so_m = m * (m - 1) // 2
print(f"The dimension of SO({m}) is calculated as: {m} * ({m} - 1) / 2 = {dim_so_m}")

# Step 3: Calculate the dimension of the unitary group U(k).
# dim(U(k)) = k^2
dim_u_k = k**2
print(f"The dimension of U({k}) is calculated as: {k}^2 = {dim_u_k}")
print("-" * 20)

# Step 4: Calculate the dimension of one component space SO(m)/U(k).
dim_component = dim_so_m - dim_u_k
print(f"The dimension of one component space SO({m})/U({k}) is: {dim_so_m} - {dim_u_k} = {dim_component}")

# Step 5: The total number of variables is the dimension of the full manifold, which has two identical components.
total_dim = 2 * dim_component
print(f"The total number of non-Grassmann variables is twice this value.")
print(f"Total variables = 2 * (dim(SO({m})) - dim(U({k}))) = 2 * ({dim_so_m} - {dim_u_k}) = 2 * {dim_component} = {total_dim}")
print("-" * 20)
print(f"The final number of non-Grassmann variables is {total_dim}.")
