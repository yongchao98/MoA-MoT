import math

def main():
    """
    Calculates the number of non-Grassmann variables for the supersymmetric
    sigma-model in symmetry class D with two replicas.
    """
    # Number of replicas as specified in the problem
    num_replicas = 2

    print("The number of non-Grassmann (bosonic) variables for symmetry class D is")
    print("given by the dimension of the symmetric space O(2N) / (O(N) x O(N)),")
    print("where N is the number of replicas.")
    print("-" * 50)

    # Set N
    N = num_replicas
    print(f"For this problem, the number of replicas N = {N}.")

    # Define the formula for the dimension of the orthogonal group O(n)
    def dim_O(n):
        return n * (n - 1) // 2

    # Calculate the dimension of the total group G = O(2N)
    n_G = 2 * N
    dim_G = dim_O(n_G)
    print(f"\nThe dimension of the total group G = O(2*N) = O({n_G}) is:")
    print(f"dim(G) = {n_G} * ({n_G} - 1) / 2 = {dim_G}")

    # Calculate the dimension of the subgroup H = O(N) x O(N)
    n_H = N
    dim_H_single = dim_O(n_H)
    dim_H = dim_H_single + dim_H_single
    print(f"\nThe dimension of the subgroup H = O(N) x O(N) = O({n_H}) x O({n_H}) is:")
    print(f"dim(H) = dim(O({n_H})) + dim(O({n_H})) = {dim_H_single} + {dim_H_single} = {dim_H}")

    # Calculate the final dimension of the manifold
    num_variables = dim_G - dim_H
    print("\nThe number of variables is dim(G) - dim(H).")
    print(f"\nFinal Equation: {dim_G} - {dim_H} = {num_variables}")
    print(f"\nTherefore, {num_variables} non-Grassmann variables are needed.")


if __name__ == "__main__":
    main()