def solve_crystal_problem():
    """
    This function calculates the dimension of the fiber and the number of energy coefficients
    for the described crystal model.

    Part A: Dimension of the fibers of pi
    The field C is identified as the 3x3 distortion tensor C_ij.
    The fiber of the bundle pi represents all possible values of C at a point.
    A 3x3 matrix has 9 independent components.
    """
    dim_fiber = 3 * 3
    
    """
    Part B: Number of coefficients specifying E
    The energy E must detect dislocations, so it depends on the dislocation density tensor α ~ curl(C).
    E must be a non-negative homogeneous polynomial of the lowest degree, implying E is quadratic in α.
    E = K_ijkl * α_ij * α_kl.
    The number of coefficients is the number of independent components of the tensor K
    for a material with cubic symmetry. This is equivalent to counting the number of
    quadratic invariants of α under the cubic group O_h.
    The 9-dimensional space of the α tensor decomposes into 4 irreducible representations for the cubic group.
    Each irreducible representation corresponds to one independent coefficient in the energy functional.
    """
    num_coeffs = 4

    # The problem asks for the output in the format "A_answer B_answer"
    # where A_answer is the dimension and B_answer is the number of coefficients.
    # The final equation is effectively the pair of answers.
    print(f"{dim_fiber} {num_coeffs}")

solve_crystal_problem()