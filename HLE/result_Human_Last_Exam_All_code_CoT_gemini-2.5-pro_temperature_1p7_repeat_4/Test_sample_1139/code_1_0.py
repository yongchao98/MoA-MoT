import math

def calculate_bosonic_dimension():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model of symmetry class D with two replicas.
    """
    # Number of replicas
    N = 2
    print(f"The number of replicas is N = {N}.")

    # The target space is M = G/H = OSp(4N|4N) / (OSp(2N|2N) x OSp(2N|2N))
    # The number of bosonic variables is dim_bos(G) - dim_bos(H)

    # Formula for bosonic dimension of OSp(2p|2q) is p*(2*p - 1) + q*(2*q + 1)
    
    # 1. Calculate the dimension of the numerator group G = OSp(4N|4N)
    p_G = 2 * N
    q_G = 2 * N
    # For N=2, G = OSp(8|8), so p_G = 4, q_G = 4
    dim_so_G = p_G * (2 * p_G - 1)
    dim_sp_G = q_G * (2 * q_G + 1)
    dim_bos_G = dim_so_G + dim_sp_G
    
    print("\nCalculating the bosonic dimension of the numerator group G = OSp(4N|4N) = OSp(8|8):")
    print(f"  p = {p_G}, q = {q_G}")
    print(f"  dim(SO(2p)) = p(2p-1) = {p_G} * (2*{p_G} - 1) = {dim_so_G}")
    print(f"  dim(Sp(2q)) = q(2q+1) = {q_G} * (2*{q_G} + 1) = {dim_sp_G}")
    print(f"  dim_bos(G) = {dim_so_G} + {dim_sp_G} = {dim_bos_G}")

    # 2. Calculate the dimension of the denominator subgroup H = OSp(2N|2N) x OSp(2N|2N)
    p_H_factor = N
    q_H_factor = N
    # For N=2, a factor of H is OSp(4|4), so p_H = 2, q_H = 2
    dim_so_H_factor = p_H_factor * (2 * p_H_factor - 1)
    dim_sp_H_factor = q_H_factor * (2 * q_H_factor + 1)
    dim_bos_H_factor = dim_so_H_factor + dim_sp_H_factor
    dim_bos_H = 2 * dim_bos_H_factor

    print("\nCalculating the bosonic dimension of the denominator group H = OSp(2N|2N) x OSp(2N|2N) = OSp(4|4) x OSp(4|4):")
    print("  First, for one factor OSp(4|4):")
    print(f"    p = {p_H_factor}, q = {q_H_factor}")
    print(f"    dim(SO(2p)) = {p_H_factor} * (2*{p_H_factor} - 1) = {dim_so_H_factor}")
    print(f"    dim(Sp(2q)) = {q_H_factor} * (2*{q_H_factor} + 1) = {dim_sp_H_factor}")
    print(f"    dim_bos(OSp(4|4)) = {dim_so_H_factor} + {dim_sp_H_factor} = {dim_bos_H_factor}")
    print(f"  The total dimension for H is twice this value:")
    print(f"  dim_bos(H) = 2 * {dim_bos_H_factor} = {dim_bos_H}")
    
    # 3. Final calculation
    num_variables = dim_bos_G - dim_bos_H
    
    print("\nThe number of non-Grassman variables is the difference between these dimensions:")
    print(f"Final Result = dim_bos(G) - dim_bos(H)")
    print(f"Final Result = {dim_bos_G} - {dim_bos_H} = {num_variables}")

calculate_bosonic_dimension()