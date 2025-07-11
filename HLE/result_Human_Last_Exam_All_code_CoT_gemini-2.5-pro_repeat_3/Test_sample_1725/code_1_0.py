def solve_crystal_model_parameters():
    """
    Calculates the two requested parameters based on the physics and symmetry of the crystal model.
    A. Dimension of the fibers of the bundle π.
    B. Number of coefficients specifying the energy functional E.
    """

    # Part A: Fiber Dimension
    # The z-axis translational invariance constrains the connection to C = A(x, y)dz.
    # The fiber is the space of possible values for the 3x3 matrix A.
    # The dimension of the space of 3x3 real matrices is 3 * 3.
    dim_fiber = 3 * 3

    # Part B: Number of Energy Coefficients
    # The energy E is a quadratic form in the derivatives of A, respecting cubic symmetry.
    # This reduces to counting the number of quadratic invariants of a 3x3 matrix
    # under the cubic group O_h. This is equivalent to the number of irreducible
    # representations into which the space of 3x3 matrices decomposes under O_h.
    # The decomposition is A_1g + T_1g + E_g + T_2g.
    # Each of the 4 irreps provides one quadratic invariant.
    num_coeffs = 4

    # Print the final result in the specified format.
    print(f"{dim_fiber} {num_coeffs}")

solve_crystal_model_parameters()