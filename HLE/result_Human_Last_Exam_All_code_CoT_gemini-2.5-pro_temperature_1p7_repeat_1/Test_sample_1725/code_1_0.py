def solve_crystal_physics_problem():
    """
    Calculates the two values requested by the user based on physics principles.

    A. Dimension of the fibers of the vector bundle pi.
    B. Number of coefficients specifying the energy functional E.
    """

    # Part A: Dimension of the fiber of the bundle pi
    # The connection C is the distortion tensor, which is a general 3x3 matrix.
    # The number of components in a 3x3 matrix is 3 * 3.
    dimension_A = 3 * 3

    # Part B: Number of coefficients for the energy density E
    # The energy E is a quadratic form in the dislocation density tensor alpha.
    # We need to find the number of independent invariants of a general
    # second-rank tensor under cubic symmetry. This can be found by
    # decomposing the tensor into its symmetric and anti-symmetric parts.

    # For the symmetric part, the number of independent coefficients
    # is the same as the number of elastic constants for a cubic crystal.
    coeffs_symmetric_part = 3

    # For the anti-symmetric part (equivalent to a pseudovector), the quadratic
    # form must be isotropic under cubic symmetry, giving 1 coefficient.
    coeffs_antisymmetric_part = 1

    # The total number of coefficients is the sum of the two parts.
    coeffs_B = coeffs_symmetric_part + coeffs_antisymmetric_part

    # The final answer requires printing both numbers separated by a space.
    print(f"{dimension_A} {coeffs_B}")

# Execute the function to display the answer.
solve_crystal_physics_problem()