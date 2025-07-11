import math

# Number of replicas
N = 2

# This script calculates the number of non-Grassman (bosonic) variables
# for the supersymmetric sigma-model of disordered systems in symmetry class D
# with N=2 replicas. The number corresponds to the dimension of the bosonic
# sector of the super-coset space G/H, where G=Osp(2N, 2N) and
# H=Osp(N, N) x Osp(N, N).

def dim_bosonic_osp(M, K):
    """
    Calculates the dimension of the bosonic part of Osp(M, 2K).
    dim_B(Osp(M, 2K)) = dim(O(M)) + dim(Sp(2K))
    """
    dim_O = M * (M - 1) / 2
    dim_Sp = K * (2 * K + 1)
    return dim_O + dim_Sp, dim_O, dim_Sp

# Step 1: Calculate the bosonic dimension of G = Osp(2N, 2N) = Osp(4, 4)
M_G = 2 * N
K_G = N
dim_B_G, dim_O_G, dim_Sp_G = dim_bosonic_osp(M_G, K_G)

# Step 2: Calculate the bosonic dimension for one component of H = Osp(N, N) = Osp(2, 2)
M_H_single = N
K_H_single = N / 2
dim_B_H_single, dim_O_H_single, dim_Sp_H_single = dim_bosonic_osp(M_H_single, K_H_single)

# The total dimension of H is twice the dimension of a single component
dim_B_H = 2 * dim_B_H_single

# Step 3: Calculate the final dimension of the bosonic coset space
result = dim_B_G - dim_B_H

# --- Output the results ---
print(f"The number of bosonic variables for symmetry class D with N={N} replicas is given by the dimension of the bosonic part of the coset space Osp({2*N}, {2*N}) / (Osp({N}, {N}) x Osp({N}, {N})).")
print(f"This dimension is calculated as: dim_B(Osp({2*N}, {2*N})) - 2 * dim_B(Osp({N}, {N})).")
print("-" * 30)

print(f"1. Dimension of the bosonic part of G = Osp({M_G}, {2*K_G}):")
print(f"   dim_B(G) = dim(O({M_G})) + dim(Sp({2*K_G}))")
print(f"            = {M_G}*({M_G}-1)/2 + {K_G}*(2*{K_G}+1)")
print(f"            = {int(dim_O_G)} + {int(dim_Sp_G)} = {int(dim_B_G)}")
print()

print(f"2. Dimension of the bosonic part of H = Osp({M_H_single}, {int(2*K_H_single)}) x Osp({M_H_single}, {int(2*K_H_single)}):")
print(f"   For a single Osp({M_H_single}, {int(2*K_H_single)}) component:")
print(f"     dim_B = dim(O({M_H_single})) + dim(Sp({int(2*K_H_single)}))")
print(f"             = {M_H_single}*({M_H_single}-1)/2 + {int(K_H_single)}*(2*{int(K_H_single)}+1)")
print(f"             = {int(dim_O_H_single)} + {int(dim_Sp_H_single)} = {int(dim_B_H_single)}")
print(f"   Total dim_B(H) = 2 * {int(dim_B_H_single)} = {int(dim_B_H)}")
print("-" * 30)

print("Final Calculation:")
print(f"Number of non-Grassman variables = dim_B(G) - dim_B(H)")
print(f"                               = {int(dim_B_G)} - {int(dim_B_H)} = {int(result)}")