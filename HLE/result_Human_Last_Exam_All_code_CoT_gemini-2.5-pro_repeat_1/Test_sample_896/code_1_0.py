import sympy

def solve_knot_polynomial_difference():
    """
    This function calculates the difference in the z^2 coefficients of the Conway polynomials
    for the closure of a given braid beta and the knot 10_4.
    """
    # Define the symbolic variable t
    t = sympy.Symbol('t')
    n = 5  # Braid group B_5

    # The braid word for beta
    braid_word_str = "s4i s4i s3i s4 s3i s2 s1i s3i s2i s2i s2i s1i"
    
    # Define the unreduced Burau representation matrices for B_5 generators
    mats = {}
    for i in range(1, n):  # i from 1 to 4
        # Matrix for sigma_i
        Si = sympy.eye(n)
        Si[i - 1, i - 1] = 1 - t
        Si[i - 1, i] = t
        Si[i, i - 1] = 1
        Si[i, i] = 0
        mats[f"s{i}"] = Si

        # Matrix for sigma_i inverse
        Sii = sympy.eye(n)
        Sii[i - 1, i - 1] = 0
        Sii[i - 1, i] = 1
        Sii[i, i - 1] = 1 / t
        Sii[i, i] = 1 - 1 / t
        mats[f"s{i}i"] = Sii

    # Compute the braid matrix B_5(beta) by multiplying generator matrices
    braid_word_list = braid_word_str.split()
    braid_matrix = sympy.eye(n)
    for gen_str in braid_word_list:
        braid_matrix = braid_matrix * mats[gen_str]

    # Calculate the non-symmetrized Alexander polynomial for the closure of beta
    # Formula: Delta(t) = det(I - B(beta)) / (1-t)
    identity_matrix = sympy.eye(n)
    determinant_val = (identity_matrix - braid_matrix).det()
    # simplify() is crucial to cancel the (1-t) factor from the determinant
    alex_poly_non_sym = sympy.simplify(determinant_val / (1 - t))

    # Symmetrize the Alexander polynomial.
    # A Laurent polynomial P(t) is symmetrized by multiplying by t^k such that P(t) = P(t^-1).
    # First, make it a polynomial by clearing the denominator.
    numer, denom = sympy.fraction(sympy.expand(alex_poly_non_sym))
    p_num = sympy.poly(numer, t)
    p_den = sympy.poly(denom, t)
    
    # Get degrees to find the normalization factor
    low_deg = p_num.degree() - p_den.degree()
    high_deg = min(p_num.monoms())[0] - max(p_den.monoms())[0]

    sym_power = -(low_deg + high_deg) / 2
    alex_poly_sym = sympy.expand(alex_poly_non_sym * (t**sym_power))
    
    # Extract coefficients a_k for the symmetric form Delta(t) = a_0 + sum(a_k * (t^k + t^-k))
    # For k > 0, a_k is the coefficient of t^k in the symmetric polynomial.
    coeffs_dict = alex_poly_sym.as_coefficients_dict()
    max_deg = 0
    for term in coeffs_dict.keys():
        if term.is_Pow and term.base == t:
            max_deg = max(max_deg, term.exp)

    coeffs_ak_beta = {}
    for k in range(1, int(max_deg) + 1):
        coeffs_ak_beta[k] = alex_poly_sym.coeff(t, k)
    
    # Calculate the z^2 coefficient of the Conway polynomial for beta-bar
    # Formula: Coeff(z^2) = sum_{k>0} k^2 * a_k
    z2_coeff_beta = sum(k*k * ak for k, ak in coeffs_ak_beta.items())

    # Now for the knot 10_4
    # The Alexander polynomial for 10_4 is Delta(t) = 2*t^2 - 5*t + 7 - 5*t^-1 + 2*t^-2.
    # This can be written as 7 - 5*(t + t^-1) + 2*(t^2 + t^-2).
    # So, a_1 = -5 and a_2 = 2.
    coeffs_ak_10_4 = {1: -5, 2: 2}
    z2_coeff_10_4 = (1**2 * coeffs_ak_10_4[1]) + (2**2 * coeffs_ak_10_4[2])
    
    # Calculate the final difference
    difference = z2_coeff_beta - z2_coeff_10_4

    # Print the results step-by-step
    print(f"The Alexander polynomial for the closure of beta, after symmetrization, is:")
    print(f"Delta_beta(t) = {sympy.pretty(alex_poly_sym)}")
    
    calc_str_beta = " + ".join([f"({k}^2)*({ak})" for k, ak in coeffs_ak_beta.items() if ak != 0])
    print(f"\nThe z^2 coefficient of the Conway polynomial for the closure of beta is:")
    print(f"c_beta = {calc_str_beta} = {z2_coeff_beta}")

    print("\nFor the knot 10_4:")
    calc_str_10_4 = f"(1^2)*({coeffs_ak_10_4[1]}) + (2^2)*({coeffs_ak_10_4[2]})"
    print(f"The z^2 coefficient of the Conway polynomial is:")
    print(f"c_10_4 = {calc_str_10_4} = {coeffs_ak_10_4[1]} + {4*coeffs_ak_10_4[2]} = {z2_coeff_10_4}")

    print(f"\nThe difference in the z^2 coefficients is:")
    print(f"c_beta - c_10_4 = {z2_coeff_beta} - {z2_coeff_10_4} = {difference}")
    print(f"<<<{difference}>>>")

solve_knot_polynomial_difference()