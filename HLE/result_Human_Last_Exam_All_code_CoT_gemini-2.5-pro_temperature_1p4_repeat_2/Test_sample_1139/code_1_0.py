def calculate_d_class_variables():
    """
    Calculates the number of non-Grassman variables for the supersymmetric
    sigma-model in symmetry class D with a given number of replicas.
    """
    # Number of replicas as specified by the user
    n = 2
    print(f"Symmetry Class: D")
    print(f"Number of replicas (n): {n}\n")

    print("The number of non-Grassman (bosonic) variables is the dimension of the bosonic sector of the target space.")
    print("This dimension is the sum of contributions from the boson-boson (BB) and fermion-fermion (FF) blocks.\n")

    # Contribution from the boson-boson block manifold: O(n,n) / (O(n) x O(n))
    dim_bb = n**2
    print(f"The dimension of the BB sector is n^2.")
    print(f"dim_BB = {n}^2 = {dim_bb}")

    # Contribution from the fermion-fermion block manifold: Sp(2n, R) / U(n)
    dim_ff = n * (n + 1)
    print(f"\nThe dimension of the FF sector is n(n+1).")
    print(f"dim_FF = {n} * ({n} + 1) = {dim_ff}")

    # Total number of bosonic variables
    total_dim = dim_bb + dim_ff
    print("\nThe total number of non-Grassman variables is the sum of these two dimensions.")
    print(f"Total Dimension = dim_BB + dim_FF")
    print(f"Result: {dim_bb} + {dim_ff} = {total_dim}")


if __name__ == "__main__":
    calculate_d_class_variables()
<<<10>>>