import math

def dim_o(p):
    """Calculates the dimension of the orthogonal Lie algebra o(p)."""
    return p * (p - 1) // 2

def dim_sp(n_2):
    """Calculates the dimension of the symplectic Lie algebra sp(2n), where n_2 = 2n."""
    if n_2 % 2 != 0:
        raise ValueError("Argument for sp must be an even number.")
    n = n_2 // 2
    return n * (n_2 + 1)

def dim_bosonic_osp(p, n_2):
    """Calculates the dimension of the bosonic part of the OSp(p|2n) Lie superalgebra."""
    dim_ortho = dim_o(p)
    dim_sympl = dim_sp(n_2)
    return dim_ortho + dim_sympl

# Number of replicas
N = 2

# For the supergroup G = OSp(2N|2N) = OSp(4|4)
p_G = 2 * N
n_2_G = 2 * N
dim_o_G = dim_o(p_G)
dim_sp_G = dim_sp(n_2_G)
dim_B_G = dim_bosonic_osp(p_G, n_2_G)

print(f"Symmetry class D with {N} replicas corresponds to the symmetric superspace OSp({2*N}|{2*N}) / (OSp({N}|{N}) x OSp({N}|{N})).")
print("\nStep 1: Calculate the bosonic dimension of G = OSp(4|4).")
print(f"The bosonic subgroup of G is O({p_G}) x Sp({n_2_G}).")
print(f"Dimension of o({p_G}) = {p_G}*({p_G}-1)/2 = {dim_o_G}")
print(f"Dimension of sp({n_2_G}) = ({n_2_G}//2)*( {n_2_G}+1) = {dim_sp_G}")
print(f"Total bosonic dimension of G, dim_B(G) = {dim_o_G} + {dim_sp_G} = {dim_B_G}")

# For the subgroup K = OSp(N|N) x OSp(N|N) = OSp(2|2) x OSp(2|2)
p_K = N
n_2_K = N
dim_o_K_single = dim_o(p_K)
dim_sp_K_single = dim_sp(n_2_K)
dim_B_K_single = dim_bosonic_osp(p_K, n_2_K)
dim_B_K = 2 * dim_B_K_single

print("\nStep 2: Calculate the bosonic dimension of K = OSp(2|2) x OSp(2|2).")
print(f"The bosonic subgroup of a single OSp({p_K}|{n_2_K}) factor is O({p_K}) x Sp({n_2_K}).")
print(f"Dimension of o({p_K}) = {p_K}*({p_K}-1)/2 = {dim_o_K_single}")
print(f"Dimension of sp({n_2_K}) = ({n_2_K}//2)*({n_2_K}+1) = {dim_sp_K_single}")
print(f"Bosonic dimension of one OSp(2|2) factor = {dim_o_K_single} + {dim_sp_K_single} = {dim_B_K_single}")
print(f"Total bosonic dimension of K, dim_B(K) = 2 * {dim_B_K_single} = {dim_B_K}")

# Final Result
result = dim_B_G - dim_B_K
print("\nStep 3: Calculate the number of non-Grassman variables.")
print("This is the difference between the bosonic dimensions of G and K.")
print(f"Number of variables = dim_B(G) - dim_B(K) = {dim_B_G} - {dim_B_K} = {result}")
