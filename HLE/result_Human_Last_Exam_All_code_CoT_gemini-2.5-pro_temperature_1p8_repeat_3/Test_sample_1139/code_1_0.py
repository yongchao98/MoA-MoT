import math

def calculate_variables_class_d(n):
    """
    Calculates the number of non-Grassmann variables for the bosonic sector of the
    supersymmetric sigma-model for disordered systems of symmetry class D with n replicas.

    This corresponds to the dimension of the symmetric space O(2n) / (O(n) x O(n)).
    """
    # Formula for the dimension of the orthogonal group O(N)
    def dim_O(N):
        return N * (N - 1) // 2

    # Step 1: Calculate the dimension of the full group G = O(2n)
    N_G = 2 * n
    dim_G = dim_O(N_G)

    # Step 2: Calculate the dimension of the subgroup K = O(n) x O(n)
    N_K_subgroup = n
    # Dimension of K is dim(O(n)) + dim(O(n))
    dim_K = 2 * dim_O(N_K_subgroup)

    # Step 3: The number of variables is the dimension of the manifold, dim(G) - dim(K)
    num_variables = dim_G - dim_K

    # --- Outputting the results as requested ---
    print("This problem asks for the number of bosonic variables for a sigma-model.")
    print("For symmetry class D with n replicas, the target space is the symmetric space O(2n)/(O(n) x O(n)).")
    print("The number of variables is the dimension of this space: dim(O(2n)) - (dim(O(n)) + dim(O(n))).")
    print("\n--- Calculation Steps ---")
    print(f"Given number of replicas, n = {n}")
    
    print(f"\n1. Dimension of the group G = O(2n) = O({N_G}):")
    print(f"   dim(O({N_G})) = {N_G} * ({N_G} - 1) / 2 = {dim_G}")
    
    print(f"\n2. Dimension of the subgroup K = O(n) x O(n) = O({N_K_subgroup}) x O({N_K_subgroup}):")
    dim_O_n = dim_O(N_K_subgroup)
    print(f"   First, dim(O({N_K_subgroup})) = {N_K_subgroup} * ({N_K_subgroup} - 1) / 2 = {dim_O_n}")
    print(f"   Then, dim(K) = dim(O({N_K_subgroup})) + dim(O({N_K_subgroup})) = {dim_O_n} + {dim_O_n} = {dim_K}")

    print("\n3. The total number of non-Grassman variables is dim(G) - dim(K).")
    print("\nThe final equation is:")
    print(f"{dim_G} - {dim_K} = {num_variables}")
    
    return num_variables

# Number of replicas given in the problem
num_replicas = 2
result = calculate_variables_class_d(num_replicas)
