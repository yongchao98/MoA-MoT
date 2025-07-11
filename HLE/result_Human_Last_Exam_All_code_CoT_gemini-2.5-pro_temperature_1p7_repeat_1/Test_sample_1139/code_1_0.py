import math

# Number of replicas as specified in the problem
n = 2

# Step 1: Calculate the dimension of the bosonic subgroup G_b = Sp(4n, R) x SO(4n)

# For Sp(4n, R), the parameter in the formula dim(Sp(2p,R)) is p = 2n
p_g = 2 * n
# For SO(4n), the parameter in the formula dim(SO(m)) is m = 4n
m_g = 4 * n

# Calculate dimension of each part
dim_sp_g = p_g * (2 * p_g + 1)
dim_so_g = m_g * (m_g - 1) // 2

# Total dimension of G_b
dim_g_b = dim_sp_g + dim_so_g

print(f"For n={n}:")
print(f"The group G = OSp({2*m_g}|{2*m_g})")
print(f"The bosonic subgroup G_b = Sp({2*p_g}, R) x SO({m_g})")
print(f"Dimension of G_b = dim(Sp({2*p_g}, R)) + dim(SO({m_g})) = {dim_sp_g} + {dim_so_g} = {dim_g_b}")
print("-" * 20)

# Step 2: Calculate the dimension of the bosonic subgroup K_b = [Sp(2n, R) x SO(2n)] x [Sp(2n, R) x SO(2n)]

# First, find the dimension of one factor, Sp(2n, R) x SO(2n)
# For Sp(2n, R), p = n
p_k = n
# For SO(2n), m = 2n
m_k = 2 * n

# Calculate dimension of each part of the factor
dim_sp_k = p_k * (2 * p_k + 1)
dim_so_k = m_k * (m_k - 1) // 2

# Dimension of one factor of K_b
dim_k_b_factor = dim_sp_k + dim_so_k

# Total dimension of K_b is twice the dimension of one factor
dim_k_b = 2 * dim_k_b_factor

print(f"The subgroup K = OSp({2*m_k}|{2*m_k}) x OSp({2*m_k}|{2*m_k})")
print(f"The bosonic subgroup K_b = [Sp({2*p_k}, R) x SO({m_k})] x [Sp({2*p_k}, R) x SO({m_k})]")
print(f"Dimension of one factor of K_b = dim(Sp({2*p_k}, R)) + dim(SO({m_k})) = {dim_sp_k} + {dim_so_k} = {dim_k_b_factor}")
print(f"Total dimension of K_b = 2 * {dim_k_b_factor} = {dim_k_b}")
print("-" * 20)

# Step 3: Calculate the final number of non-Grassman variables
num_variables = dim_g_b - dim_k_b

print("Final Calculation:")
print(f"Number of non-Grassman variables = dim(G_b) - dim(K_b)")
print(f"The number of variables is {dim_g_b} - {dim_k_b} = {num_variables}")