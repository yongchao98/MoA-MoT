def solve_crystal_problem():
    """
    This function calculates the answers to the two parts of the problem based on
    the physical and mathematical interpretation of the continuum crystal model.
    """

    # Part A: Dimension of the fiber of the bundle pi
    # The field C is interpreted as the 3x3 distortion tensor. The fiber of the
    # bundle is the space of possible values of C at a point, which is the
    # space of 3x3 matrices.
    # The dimension of the space of 3x3 matrices is 3 * 3.
    dimension_of_fiber = 3 * 3

    # Part B: Number of coefficients specifying the energy E
    # The energy E is a quadratic form in the dislocation density tensor alpha,
    # which is a rank-2 tensor (a 3x3 matrix). The number of coefficients is
    # the number of independent quadratic invariants of a rank-2 tensor under
    # cubic symmetry (O_h group). This equals the number of irreducible
    # representations in the decomposition of the space of rank-2 tensors.
    # The decomposition is A_1g + T_1g + E_g + T_2g.
    # There are 4 irreducible representations.
    number_of_coefficients = 4

    # The problem asks for the output in the format "A B"
    print(f"{dimension_of_fiber} {number_of_coefficients}")

solve_crystal_problem()