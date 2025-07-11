def calculate_bosonic_variables():
    """
    Calculates the number of non-Grassman variables for the bosonic sector
    of the supersymmetric sigma-model for symmetry class D.
    """
    # The number of replicas
    N = 2

    # For symmetry class D, the dimensionality of the bosonic sector of the
    # target supermanifold is given by the formula 2 * N^2.
    num_variables = 2 * N**2

    # Print the explanation and the final equation with all numbers.
    print("The number of non-Grassman variables needed to parametrize the bosonic sector of the")
    print("supersymmetric sigma-model for disordered systems of symmetry class D with N replicas is given by the formula: 2 * N^2.")
    print(f"\nFor N = {N} replicas, the calculation is:")
    print(f"2 * {N}^2 = 2 * {N*N} = {num_variables}")

if __name__ == "__main__":
    calculate_bosonic_variables()