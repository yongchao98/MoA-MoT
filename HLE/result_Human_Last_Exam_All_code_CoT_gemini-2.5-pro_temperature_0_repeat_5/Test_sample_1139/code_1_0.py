import math

def main():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model of symmetry class D with a given number of replicas.
    """
    # Number of replicas as specified in the problem
    n = 2

    print("To find the number of non-Grassman variables for the supersymmetric sigma-model of symmetry class D, we must calculate the dimension of its bosonic target manifold.")
    print(f"For class D with n={n} replicas, this manifold is the symmetric space G/K = O(2n) / (O(n) x O(n)).")
    print("The number of variables is the dimension of this space, which is calculated as: dim(G) - dim(K).")
    print("-" * 40)

    # --- Calculation ---
    # The dimension of the orthogonal group O(N) is N*(N-1)/2.

    # Step 1: Calculate the dimension of G = O(2n)
    N_G = 2 * n
    dim_G = N_G * (N_G - 1) // 2

    print(f"Step 1: Calculate the dimension of the group G = O(2n) = O({N_G}).")
    print(f"dim(O({N_G})) = ({N_G} * ({N_G} - 1)) / 2 = {dim_G}")
    print("-" * 40)

    # Step 2: Calculate the dimension of the subgroup K = O(n) x O(n)
    N_K = n
    dim_O_n = N_K * (N_K - 1) // 2
    dim_K = dim_O_n + dim_O_n

    print(f"Step 2: Calculate the dimension of the subgroup K = O(n) x O(n) = O({N_K}) x O({N_K}).")
    print(f"The dimension of a single O({N_K}) group is ({N_K} * ({N_K} - 1)) / 2 = {dim_O_n}.")
    print(f"dim(K) = dim(O({N_K})) + dim(O({N_K})) = {dim_O_n} + {dim_O_n} = {dim_K}")
    print("-" * 40)

    # Step 3: Calculate the final number of variables
    num_variables = dim_G - dim_K

    print("Step 3: Calculate the final number of variables by finding the dimension of the manifold G/K.")
    print(f"Number of variables = dim(G) - dim(K)")
    print(f"The final equation is: {dim_G} - {dim_K} = {num_variables}")
    print("-" * 40)
    print(f"The number of non-Grassman variables needed is {num_variables}.")

if __name__ == "__main__":
    main()