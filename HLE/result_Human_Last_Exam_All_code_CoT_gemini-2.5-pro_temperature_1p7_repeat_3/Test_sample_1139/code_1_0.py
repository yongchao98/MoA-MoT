import math

def solve_class_d_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric sigma-model
    of symmetry class D with a given number of replicas.
    """
    # Number of replicas as specified in the problem
    N = 2

    # --- Step 1: Calculate the dimension of the bosonic part of G = Osp(2N|2N) ---
    # The bosonic subgroup of Osp(2N|2N) is SO(2N) x Sp(2N).
    # Its dimension is dim(SO(2N)) + dim(Sp(2N)).

    # Dimension of SO(k) = k * (k - 1) / 2
    dim_so_2N = (2 * N) * (2 * N - 1) // 2

    # Dimension of Sp(2k) = k * (2k + 1)
    # Here, the group is Sp(2N), so the parameter for the formula is k=N.
    dim_sp_2N = N * (2 * N + 1)

    # Total bosonic dimension of G = Osp(2N|2N)
    dim_G_bos = dim_so_2N + dim_sp_2N

    # --- Step 2: Calculate the dimension of the bosonic part of H = U(N|N) ---
    # The bosonic subgroup of U(N|N) is U(N) x U(N).
    # Its dimension is dim(U(N)) + dim(U(N)).

    # Dimension of U(k) = k^2
    dim_u_N = N * N

    # Total bosonic dimension of H = U(N|N)
    dim_H_bos = dim_u_N + dim_u_N

    # --- Step 3: Calculate the final number of variables ---
    # This is the dimension of the bosonic part of the coset space G/H.
    num_bosonic_variables = dim_G_bos - dim_H_bos

    # --- Step 4: Print the detailed calculation ---
    print("This script calculates the number of non-Grassman variables for the bosonic sector of the supersymmetric sigma-model for disordered systems of symmetry class D with N replicas.")
    print(f"The number of replicas is N = {N}.\n")
    print("The target superspace for this model is G/H = Osp(2N|2N) / U(N|N).")
    print("The number of bosonic variables is dim_bos(G) - dim_bos(H).\n")
    
    print("--- Calculation ---")
    print(f"1. Bosonic dimension of G = Osp({2*N}|{2*N}):")
    print(f"   dim_bos(G) = dim(SO({2*N})) + dim(Sp({2*N}))")
    print(f"   dim_bos(G) = {dim_so_2N} + {dim_sp_2N} = {dim_G_bos}\n")

    print(f"2. Bosonic dimension of H = U({N}|{N}):")
    print(f"   dim_bos(H) = dim(U({N})) + dim(U({N}))")
    print(f"   dim_bos(H) = {dim_u_N} + {dim_u_N} = {dim_H_bos}\n")

    print("3. Number of non-Grassman (bosonic) variables:")
    print("   Number = dim_bos(G) - dim_bos(H)")
    print(f"   Number = {dim_G_bos} - {dim_H_bos} = {num_bosonic_variables}")
    
# Run the calculation
solve_class_d_variables()
