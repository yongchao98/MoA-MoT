import math

def solve_crystal_problem():
    """
    Solves the two-part physics problem about a continuum crystal model.
    A. What is the dimension of pi's fibers?
    B. How many coefficients specify E?
    """

    # Part A: Dimension of pi's fibers
    # The connection C on TR^3 corresponds to the distortion tensor C_ij, a 3x3 matrix.
    # The fiber of the bundle pi at a point x is the space of all possible values of C(x),
    # which is the space of 3x3 real matrices. The dimension is 3 * 3.
    dim_pi_fiber = 3 * 3

    # Part B: Number of coefficients specifying E
    # The energy density E is an O_h-invariant quadratic form on the space of 3-tensors
    # X_kij = d_k C_ij. The number of coefficients is the dimension of the space of
    # such invariant forms, calculated using character theory.
    # N = (1 / (2*|G|)) * sum_{g in G} [ (chi_v(g)^6) + (chi_v(g^2))^3 ]
    # where v is the 3D vector representation of G = O_h.

    G_order = 48
    
    # Conjugacy classes data for O_h group: [class_size, chi_v(g), chi_v(g^2)]
    # g^2 for an element in a class C can fall into another class C'. The value
    # chi_v(g^2) is the character of the class of g^2.
    # Classes are ordered as: E, 3C2, 6C4, 8C3, 6C'2, i, 3sigma_h, 6S4, 8S6, 6sigma_d.
    classes_data = [
        [1, 3, 3],   # E -> E
        [3, -1, 3],  # C2 -> E
        [6, 1, -1],  # C4 -> C2
        [8, 0, 0],   # C3 -> C3^2 (same character as C3)
        [6, -1, 3],  # C'2 -> E
        [1, -3, 3],  # i -> E
        [3, 1, 3],   # sigma_h -> E
        [6, -1, -1], # S4 -> C2
        [8, 0, 0],   # S6 -> C3
        [6, 1, 3]    # sigma_d -> E
    ]

    total_sum_for_char_formula = 0
    for data in classes_data:
        class_size, chi_v_g, chi_v_g_squared = data
        
        # The representation on the space of 3-tensors V has character chi_V(g) = chi_v(g)^3.
        chi_V_g = chi_v_g ** 3
        chi_V_g_squared = chi_v_g_squared ** 3
        
        # The summand in the character formula for the number of quadratic invariants.
        summand = (chi_V_g ** 2) + chi_V_g_squared
        
        total_sum_for_char_formula += class_size * summand

    # The number of coefficients N is given by the formula for the dimension of the invariant subspace.
    num_coeffs = total_sum_for_char_formula / (2 * G_order)
    
    # The result must be an integer as it represents a dimension.
    num_coeffs_int = int(round(num_coeffs))

    # Print the final answer in the format "A B"
    print(f"{dim_pi_fiber} {num_coeffs_int}")

solve_crystal_problem()