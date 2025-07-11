import math

def dim_orthogonal(p):
    """Calculates the dimension of the orthogonal group O(p)."""
    return p * (p - 1) // 2

def dim_symplectic_real(p):
    """
    Calculates the dimension of the real symplectic group Sp(p,R).
    p must be an even integer. The formula is n * (2n + 1) where p=2n.
    """
    if p % 2 != 0:
        raise ValueError(f"The dimension 'p' for Sp(p,R) must be even, but got {p}.")
    n = p // 2
    return n * (2 * n + 1)

def solve_class_d_bosonic_variables():
    """
    Calculates the number of non-Grassman variables needed to parametrize the
    bosonic sector of the supersymmetric sigma-model for disordered systems of
    symmetry class D with a given number of replicas.
    """
    N = 2  # Number of replicas

    print(f"Calculating the number of bosonic variables for symmetry class D with N = {N} replicas.")
    print("-" * 70)
    print("The target space is the bosonic part of OSp(2N|2N) / (OSp(N|N) x OSp(N|N)).")
    print("This corresponds to the coset G/H where G = O(2N) x Sp(2N,R) and H = [O(N) x Sp(N,R)] x [O(N) x Sp(N,R)].")
    print("The number of variables is dim(G) - dim(H).\n")

    # Dimensions for the full group G = O(2N) x Sp(2N,R)
    group_g_dim_o = dim_orthogonal(2 * N)
    group_g_dim_sp = dim_symplectic_real(2 * N)
    total_dim_g = group_g_dim_o + group_g_dim_sp
    
    print("Step 1: Calculate the dimension of the parent group G = O(4) x Sp(4,R).")
    print(f"  - dim(O({2*N})) = {2*N}*({2*N}-1)/2 = {group_g_dim_o}")
    print(f"  - dim(Sp({2*N},R)), with n={N}, is {N}*(2*{N}+1) = {group_g_dim_sp}")
    print(f"  - Total dim(G) = dim(O({2*N})) + dim(Sp({2*N},R)) = {group_g_dim_o} + {group_g_dim_sp} = {total_dim_g}\n")

    # Dimensions for the subgroup H = [O(N)xSp(N,R)] x [O(N)xSp(N,R)]
    group_h_dim_o = dim_orthogonal(N)
    group_h_dim_sp = dim_symplectic_real(N)
    total_dim_h = 2 * (group_h_dim_o + group_h_dim_sp)

    print(f"Step 2: Calculate the dimension of the subgroup H = [O({N}) x Sp({N},R)]^2.")
    print(f"  For one component O({N}) x Sp({N},R):")
    print(f"  - dim(O({N})) = {N}*({N}-1)/2 = {group_h_dim_o}")
    n_sub = N // 2
    print(f"  - dim(Sp({N},R)), with n={n_sub}, is {n_sub}*(2*{n_sub}+1) = {group_h_dim_sp}")
    print(f"  - Total dim(H) = 2 * (dim(O({N})) + dim(Sp({N},R))) = 2 * ({group_h_dim_o} + {group_h_dim_sp}) = {total_dim_h}\n")
    
    # Final result
    num_variables = total_dim_g - total_dim_h
    print("Step 3: Calculate the number of non-Grassman variables.")
    print(f"  Number of variables = dim(G) - dim(H) = {total_dim_g} - {total_dim_h} = {num_variables}")
    
    return num_variables

if __name__ == '__main__':
    solve_class_d_bosonic_variables()
