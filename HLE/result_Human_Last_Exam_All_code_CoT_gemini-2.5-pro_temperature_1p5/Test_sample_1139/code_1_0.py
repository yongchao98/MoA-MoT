import math

def calculate_bosonic_variables():
    """
    Calculates the number of non-Grassman variables for a supersymmetric
    sigma-model in symmetry class D with a given number of replicas.
    """
    # Number of replicas, as specified in the problem
    num_replicas = 2
    
    # For symmetry class D, the bosonic sector of the sigma-model is parametrized
    # by the Riemannian symmetric space G/H, where G = O(2*N_r) and H = U(N_r),
    # with N_r being the number of replicas.

    # Step 1: Determine the specific groups G and H for N_r = 2 replicas.
    n_for_O = 2 * num_replicas
    n_for_U = num_replicas
    
    print(f"The number of bosonic variables corresponds to the dimension of the coset space O({n_for_O})/U({n_for_U}).")
    print("The dimension is calculated as: dim(O({n_for_O})) - dim(U({n_for_U}))\n")

    # Step 2: Calculate the dimension of the orthogonal group G = O(n_for_O).
    # The formula for the dimension of O(n) is n * (n - 1) / 2.
    dim_O = n_for_O * (n_for_O - 1) // 2
    print(f"First, we calculate the dimension of O({n_for_O}):")
    print(f"dim(O({n_for_O})) = ({n_for_O} * ({n_for_O} - 1)) / 2 = {dim_O}")

    # Step 3: Calculate the dimension of the unitary group H = U(n_for_U).
    # The formula for the dimension of U(n) is n^2.
    dim_U = n_for_U ** 2
    print(f"\nNext, we calculate the dimension of U({n_for_U}):")
    print(f"dim(U({n_for_U})) = {n_for_U}^2 = {dim_U}")

    # Step 4: Calculate the final dimension of the coset space G/H.
    num_variables = dim_O - dim_U
    
    print("\nFinally, we subtract the dimension of H from the dimension of G.")
    print("The final equation for the number of variables is:")
    print(f"{dim_O} - {dim_U} = {num_variables}")

if __name__ == '__main__':
    calculate_bosonic_variables()