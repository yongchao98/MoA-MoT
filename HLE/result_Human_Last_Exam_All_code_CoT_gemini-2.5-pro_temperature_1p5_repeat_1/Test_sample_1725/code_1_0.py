def solve_crystal_properties():
    """
    Calculates the dimension of the bundle's fibers and the number of coefficients for the energy density.
    """

    # Part A: Dimension of pi's fibers
    # The field C is the distortion tensor, a 3x3 matrix at each point in space.
    # The dimension of the space of 3x3 matrices is 3 * 3.
    dim_fiber = 3 * 3

    # Part B: Number of coefficients specifying E
    # The energy E is a quadratic form in the dislocation density tensor alpha.
    # E = L_ijkl * alpha_ij * alpha_kl
    # For a cubic crystal, we need to count the number of independent components
    # of the 4th-rank tensor L that are invariant under the cubic symmetry group.
    # Based on group representation theory, this number is 4.
    num_coeffs = 4

    print(f"{dim_fiber} {num_coeffs}")

solve_crystal_properties()