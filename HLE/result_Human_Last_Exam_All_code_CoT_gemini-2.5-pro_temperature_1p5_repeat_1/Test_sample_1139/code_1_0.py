import math

def main():
    """
    Calculates the number of non-Grassman variables needed to parametrize
    the bosonic sector of the supersymmetric sigma-model for disordered systems
    of symmetry class D with a given number of replicas.
    """
    # Number of replicas as given in the problem
    n = 2

    print(f"The problem is to find the number of non-Grassman (bosonic) variables for a class D SUSY sigma-model with n = {n} replicas.")
    print("This corresponds to the dimension of the bosonic part of the target supermanifold.")
    print("For class D, the target manifold is the symmetric superspace M = G/K, where:")
    print(f"G = OSp(2n|2n) = OSp({2*n}|{2*n})")
    print(f"K = OSp(n|n) x OSp(n|n) = OSp({n}|{n}) x OSp({n}|{n})")
    print("\nThe number of variables is dim_B(M) = dim_B(G) - dim_B(K).\n")

    # The bosonic dimension of OSp(2p|2q) is dim(SO(2p)) + dim(Sp(2q))
    # where dim(SO(m)) = m*(m-1)/2 and dim(Sp(2q)) = q*(2q+1)

    # --- Step 1: Calculate the bosonic dimension of G = OSp(4|4) ---
    print("--- Calculating dim_B(G) for G = OSp(4|4) ---")
    p_g = (2 * n) // 2
    q_g = (2 * n) // 2

    # Dimension of SO(2n) = SO(4)
    dim_so_g = (2 * n) * (2 * n - 1) // 2
    # Dimension of Sp(2n) = Sp(4)
    dim_sp_g = q_g * (2 * q_g + 1)
    
    dim_b_g = dim_so_g + dim_sp_g
    
    print(f"dim_B(OSp(4|4)) = dim(SO(4)) + dim(Sp(4))")
    print(f"dim_B(OSp(4|4)) = {dim_so_g} + {dim_sp_g} = {dim_b_g}\n")

    # --- Step 2: Calculate the bosonic dimension of K = OSp(2|2) x OSp(2|2) ---
    print("--- Calculating dim_B(K) for K = OSp(2|2) x OSp(2|2) ---")
    
    # We first calculate the dimension for one factor: OSp(n|n) = OSp(2|2)
    p_k = n // 2
    q_k = n // 2
    
    # Dimension of SO(n) = SO(2)
    dim_so_k_factor = n * (n - 1) // 2
    # Dimension of Sp(n) = Sp(2)
    dim_sp_k_factor = q_k * (2 * q_k + 1)
    
    dim_b_k_factor = dim_so_k_factor + dim_sp_k_factor
    
    print("For one factor OSp(2|2):")
    print(f"dim_B(OSp(2|2)) = dim(SO(2)) + dim(Sp(2))")
    print(f"dim_B(OSp(2|2)) = {dim_so_k_factor} + {dim_sp_k_factor} = {dim_b_k_factor}")

    # The total dimension of K is twice the dimension of one factor
    dim_b_k = 2 * dim_b_k_factor
    print(f"\nFor the full subgroup K = OSp(2|2) x OSp(2|2):")
    print(f"dim_B(K) = 2 * dim_B(OSp(2|2)) = 2 * {dim_b_k_factor} = {dim_b_k}\n")

    # --- Step 3: Final Calculation ---
    print("--- Final Calculation ---")
    result = dim_b_g - dim_b_k
    print("The number of non-Grassman variables is dim_B(M) = dim_B(G) - dim_B(K)")
    print(f"Number of variables = {dim_b_g} - {dim_b_k} = {result}")

if __name__ == "__main__":
    main()