import sys

def solve_class_d_bosonic_variables():
    """
    Calculates the number of non-Grassmann variables for the supersymmetric
    sigma-model in symmetry class D with two replicas.
    """
    # Number of replicas for the disordered system
    N = 2

    # The number of non-Grassmann (bosonic) variables is given by the dimension
    # of the bosonic part of the target space G/H, where G = O(2N|2N) and
    # H = O(N|N) x O(N|N). This dimension is dim_B(g) - dim_B(h).

    # Dimension of the bosonic part of the Lie superalgebra g = o(2N|2N).
    # The bosonic part is o(2N) + o(2N).
    # The dimension of o(k) is k*(k-1)/2.
    # So, dim_B(g) = 2 * (2*N * (2*N - 1) / 2) = 2*N*(2*N-1).
    dim_g_B = 2 * N * (2 * N - 1)

    # Dimension of the bosonic part of the Lie superalgebra h = o(N|N) + o(N|N).
    # The bosonic part is o(N)+o(N) + o(N)+o(N).
    # So, dim_B(h) = 4 * (N * (N - 1) / 2) = 2*N*(N-1).
    dim_h_B = 2 * N * (N - 1)

    # The number of variables is the difference.
    num_variables = dim_g_B - dim_h_B

    print(f"The number of non-Grassmann variables for the supersymmetric sigma-model in symmetry class D with {N} replicas is calculated as follows:")
    print("The target space is the supermanifold O(2N|2N) / (O(N|N) x O(N|N)).")
    print("The number of bosonic variables is the dimension of the bosonic part of this space.")
    print("This is calculated as the difference in the dimensions of the bosonic parts of the corresponding Lie superalgebras, dim_B(g) - dim_B(h).")
    print(f"\nFor N = {N}:")
    print(f"Dimension of the bosonic part of g, dim_B(g) = 2*N*(2*N-1) = 2*{N}*(2*{N}-1) = {dim_g_B}")
    print(f"Dimension of the bosonic part of h, dim_B(h) = 2*N*(N-1) = 2*{N}*({N}-1) = {dim_h_B}")
    print("\nThe final number of non-Grassmann variables is given by the equation:")
    print(f"{dim_g_B} - {dim_h_B} = {num_variables}")

solve_class_d_bosonic_variables()