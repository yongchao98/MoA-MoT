def display_generating_function_asymptotics():
    """
    Performs and displays an asymptotic analysis of the billiard generating function H(s, s').

    This function derives and prints the leading-order behavior of H(s, s')
    in the limit as the separation |s' - s| -> 0, including the effect
    of local boundary curvature κ(s).
    """

    # --- Symbolic representation of variables ---
    gen_func = "H(s, s')"
    separation = "|s' - s|"
    curvature = "κ(s)"

    # --- Introduction ---
    print("Asymptotic Analysis of the Billiard Generating Function H(s, s')")
    print("-" * 60)
    print(f"The generating function {gen_func} is the Euclidean distance between two points on the billiard boundary.")
    print(f"We are interested in its behavior as the arc-length separation, {separation}, approaches zero.\n")

    # --- Derivation of Terms ---
    print("The expansion can be written as a series:")
    print(f"{gen_func} = (Term 1) + (Term 2) + ...\n")

    # Term 1: The leading-order term
    term1_coeff = 1
    term1_expr = f"{separation}"
    print(f"Term 1: The dominant term is the arc-length separation itself.")
    print(f"  - Expression: {term1_coeff} * {term1_expr}")
    print(f"  - This represents the chord length in the limit of zero curvature (a straight line).\n")

    # Term 2: The curvature correction
    # The coefficient is -1/24, resulting from the Taylor expansion.
    # H^2 ≈ (Δs)^2 - (1/12)κ^2(Δs)^4
    # H ≈ Δs * sqrt(1 - (1/12)κ^2(Δs)^2)
    # Using sqrt(1-x) ≈ 1 - x/2, H ≈ Δs * (1 - (1/2)*(1/12)κ^2(Δs)^2) = Δs - (1/24)κ^2(Δs)^3
    term2_coeff_num = -1
    term2_coeff_den = 24
    term2_expr = f"{curvature}² * {separation}³"
    print(f"Term 2: The first-order correction due to the boundary's geometry.")
    print(f"  - This term incorporates the local curvature, {curvature}.")
    print(f"  - Coefficient: {term2_coeff_num}/{term2_coeff_den}")
    print(f"  - Expression: ({term2_coeff_num}/{term2_coeff_den}) * {term2_expr}")
    print(f"  - The negative sign indicates the chord length is shorter than the arc length for a curved boundary.\n")

    # --- Final Equation ---
    print("Combining these terms, the asymptotic expansion is:")
    print("-" * 60)
    # The final equation is built showing each number explicitly as requested.
    final_equation = (f"{gen_func} ≈ {term1_coeff} * {separation} "
                      f"+ ({term2_coeff_num} / {term2_coeff_den}) * {curvature}² * {separation}³ "
                      f"+ O({separation}⁴)")
    
    # We print a cleaner version for final presentation
    final_equation_clean = (f"{gen_func} ≈ {separation} - (1/24) * {curvature}² * {separation}³ + O({separation}⁴)")
    
    print(final_equation_clean)
    print("-" * 60)


if __name__ == "__main__":
    display_generating_function_asymptotics()