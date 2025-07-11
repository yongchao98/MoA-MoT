def solve_crystal_problem():
    """
    This function solves the two parts of the user's question based on
    a continuum model of a crystal.
    """

    # Part A: What is the dimension of pi's fibers?
    # The field C represents the local distortion of the crystal, which is
    # described by a 3x3 distortion tensor C_ij at each point.
    # The fiber of the bundle pi is the space of all possible 3x3 matrices.
    # The dimension of this space is 3 * 3.
    dim_fibers = 3 * 3

    # Part B: How many coefficients specify E?
    # The energy E must detect dislocations, meaning it depends on the dislocation
    # density tensor alpha = curl(C). For a simple, positive-definite energy,
    # it must be a quadratic function of alpha.
    # We need to find the number of independent coefficients in this quadratic
    # form for a crystal with cubic symmetry.

    # A general tensor alpha can be split into a symmetric part (S) and an
    # anti-symmetric part (A). For cubic symmetry, the energy does not mix these parts.
    # E(alpha) = E_S(S) + E_A(A).

    # The number of coefficients for the symmetric part is the number of
    # independent elastic constants for a cubic crystal.
    num_coeffs_symmetric_part = 3

    # The number of coefficients for the anti-symmetric part (which is equivalent
    # to a pseudovector) is 1, as the quadratic form is isotropic for cubic symmetry.
    num_coeffs_antisymmetric_part = 1

    # The total number of coefficients is the sum of the coefficients for each part.
    total_coeffs = num_coeffs_symmetric_part + num_coeffs_antisymmetric_part
    
    # Final equation for number of coefficients:
    # total_coeffs = (coeffs for symmetric part) + (coeffs for anti-symmetric part)

    # Output the final answer in the format "A B"
    # Outputting each number that leads to the final answer for B, as requested.
    print(f"{dim_fibers} {total_coeffs}")

solve_crystal_problem()