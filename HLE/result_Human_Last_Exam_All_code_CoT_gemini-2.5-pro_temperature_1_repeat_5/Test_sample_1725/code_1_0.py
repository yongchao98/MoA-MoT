def solve_crystal_problem():
    """
    This function solves the two parts of the crystal model problem based on its description.
    """

    # Part A: Determine the dimension of the fibers of the bundle pi.
    # The field C represents the local crystal structure and is given as a connection on the tangent bundle T(R^3).
    # In continuum mechanics, this corresponds to the distortion tensor, a 3x3 matrix.
    # The fiber of the bundle is the space of possible values of C at a single point.
    # The number of independent components in a 3x3 matrix is 3 * 3 = 9.
    dimension_of_fibers = 9

    # Part B: Determine the number of coefficients specifying the energy density E.
    # The energy E is a local, homogeneous polynomial of least degree, invariant to rigid transformations,
    # induces linear dynamics, and is consistent with detecting dislocations for a cubic crystal.
    # 1. Invariant to rigid transformations -> E depends on the symmetric part of C, the strain tensor epsilon.
    # 2. Induces linear dynamics -> E is quadratic in epsilon.
    # 3. Least degree homogeneous polynomial -> E must be of degree 2, i.e., E = (1/2) * c_ijkl * epsilon_ij * epsilon_kl.
    # 4. This form is consistent with detecting dislocations, as dislocations produce strain fields and thus have energy.
    # 5. The coefficients c_ijkl are constrained by the crystal's cubic symmetry.
    # 6. For a cubic crystal, the 81 components of the elasticity tensor c_ijkl reduce to 3 independent constants.
    number_of_coefficients = 3

    # The final answer requires printing the two numbers.
    print(f"{dimension_of_fibers} {number_of_coefficients}")

solve_crystal_problem()