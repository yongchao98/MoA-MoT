import math

def calculate_bosonic_variables():
    """
    Calculates the number of non-Grassmann variables for the supersymmetric sigma-model
    of symmetry class D with a specified number of replicas.
    """
    # Number of replicas
    n = 2

    # Formulas for Lie group dimensions
    def dim_so(k):
        """Calculates the dimension of the special orthogonal group SO(k)."""
        return k * (k - 1) // 2

    def dim_sp(k_even):
        """Calculates the dimension of the real symplectic group Sp(k_even, R)."""
        if k_even % 2 != 0:
            raise ValueError("The argument for Sp must be an even number.")
        k = k_even // 2
        return k * (2 * k + 1)

    print(f"Calculating the number of bosonic variables for symmetry class D with n = {n} replicas.")
    print("The target supermanifold is G/H, where G = OSp(2n|2n) and H = OSp(n|n) x OSp(n|n).")
    print("The number of bosonic variables equals the dimension of the bosonic base G_B/H_B.")
    print("-" * 30)

    # Calculate dimension of the bosonic supergroup G_B = SO(2n) x Sp(2n,R)
    g_so_k = 2 * n
    g_sp_k = 2 * n
    dim_g_b_so = dim_so(g_so_k)
    dim_g_b_sp = dim_sp(g_sp_k)
    dim_g_b = dim_g_b_so + dim_g_b_sp

    print(f"Dimension of the bosonic group G_B = SO({g_so_k}) x Sp({g_sp_k}, R):")
    print(f"dim(SO({g_so_k})) = {g_so_k} * ({g_so_k} - 1) / 2 = {dim_g_b_so}")
    print(f"dim(Sp({g_sp_k}, R)) = {g_sp_k//2} * (2 * {g_sp_k//2} + 1) = {dim_g_b_sp}")
    print(f"Total dim(G_B) = {dim_g_b_so} + {dim_g_b_sp} = {dim_g_b}")
    print("-" * 30)
    
    # Calculate dimension of the bosonic subgroup H_B = (SO(n) x Sp(n,R)) x (SO(n) x Sp(n,R))
    h_so_k = n
    h_sp_k = n
    dim_h_b_so_factor = dim_so(h_so_k)
    dim_h_b_sp_factor = dim_sp(h_sp_k)
    dim_h_b_factor = dim_h_b_so_factor + dim_h_b_sp_factor
    dim_h_b = 2 * dim_h_b_factor
    
    print(f"Dimension of the bosonic subgroup H_B = [SO({h_so_k}) x Sp({h_sp_k}, R)] x [SO({h_so_k}) x Sp({h_sp_k}, R)]:")
    print(f"Dimension of one [SO({h_so_k}) x Sp({h_sp_k}, R)] factor:")
    print(f"  dim(SO({h_so_k})) = {h_so_k} * ({h_so_k} - 1) / 2 = {dim_h_b_so_factor}")
    print(f"  dim(Sp({h_sp_k}, R)) = {h_sp_k//2} * (2 * {h_sp_k//2} + 1) = {dim_h_b_sp_factor}")
    print(f"Dimension of one bosonic factor of H = {dim_h_b_so_factor} + {dim_h_b_sp_factor} = {dim_h_b_factor}")
    print(f"Total dim(H_B) = 2 * {dim_h_b_factor} = {dim_h_b}")
    print("-" * 30)

    # Final result is the dimension of the coset space G_B/H_B
    num_variables = dim_g_b - dim_h_b
    
    print("The number of non-Grassmann variables is dim(G_B/H_B) = dim(G_B) - dim(H_B).")
    print(f"Result = {dim_g_b} - {dim_h_b} = {num_variables}")
    
calculate_bosonic_variables()