def solve_crystal_problem():
    """
    This function calculates and prints the answers to the two-part physics problem.

    A. What is the dimension of pi's fibers?
    The connection C is a representation of the local crystal distortion, which is
    a 3x3 matrix. The space of all 3x3 real matrices has dimension 3*3 = 9.

    B. How many coefficients specify E?
    The energy functional E must be quadratic for linear dynamics and respect cubic
    symmetry. It consists of two parts:
    1. An elastic part based on the symmetric strain tensor, which has 3
       independent coefficients for a cubic crystal.
    2. A dislocation part based on the curl of the distortion field (a general
       second-rank tensor). For cubic symmetry, there are 4 independent
       quadratic invariants that can be formed from this tensor.
    The total number of coefficients is the sum: 3 + 4 = 7.
    """

    # Part A: Dimension of pi's fibers
    dimension_A = 9

    # Part B: Number of coefficients specifying E
    elastic_coeffs = 3
    dislocation_coeffs = 4
    coefficients_B = elastic_coeffs + dislocation_coeffs
    
    # Print the answer in the required format "A B"
    print(f"{dimension_A} {coefficients_B}")

solve_crystal_problem()