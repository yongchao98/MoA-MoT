import math

def solve_sigma_model_dimension():
    """
    Calculates the number of non-Grassmann variables for the bosonic sector of the
    supersymmetric sigma-model with n=2 replicas for symmetry class D.
    """
    # Number of replicas
    n = 2

    # The bosonic sector is parametrized by the symmetric space O(2n) / (O(n) x O(n)).
    # We need to calculate the dimension of this space.

    def dim_orthogonal_group(N):
        """Calculates the dimension of the orthogonal group O(N)."""
        if N <= 0:
            return 0
        return N * (N - 1) // 2

    # Parameters for the main group G = O(2n) and subgroup K = O(n) x O(n)
    N_G = 2 * n
    N_K = n

    # Calculate the dimension of the main group G
    dim_G = dim_orthogonal_group(N_G)

    # Calculate the dimension of the subgroup K
    # dim(K) = dim(O(n)) + dim(O(n))
    dim_K_component = dim_orthogonal_group(N_K)
    dim_K_total = 2 * dim_K_component

    # The number of variables is the dimension of the quotient space G/K
    num_variables = dim_G - dim_K_total

    # Print the step-by-step calculation
    print("The number of non-Grassmann variables is the dimension of the symmetric space O(2n) / (O(n) x O(n)) for n=2.")
    print(f"This corresponds to the space O({N_G}) / (O({N_K}) x O({N_K})).")
    print("\nThe dimension is calculated as: dim(G) - dim(K)")
    print(f"dim(O({N_G})) = {N_G}*({N_G}-1)/2 = {dim_G}")
    print(f"dim(O({N_K})) = {N_K}*({N_K}-1)/2 = {dim_K_component}")
    print("\nThe final calculation is:")
    print(f"dim(O({N_G})) - (dim(O({N_K})) + dim(O({N_K}))) = {dim_G} - ({dim_K_component} + {dim_K_component}) = {num_variables}")
    print(f"\nTherefore, {num_variables} non-Grassmann variables are needed to parametrize the bosonic sector.")

solve_sigma_model_dimension()