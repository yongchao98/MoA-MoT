def generate_formula():
    """
    This function prints the derived formula for P(n).
    """

    # Coefficients for the numerator of the 1/n^2 term
    coeff_L2_n2 = 3
    coeff_L1_n2 = -2
    coeff_L0_n2 = 2
    denominator_n2 = 24

    # Coefficients for the numerator of the 1/n^3 term
    coeff_L3_n3 = 1
    coeff_L2_n3 = -2
    coeff_L1_n3 = 2
    denominator_n3 = 48

    # Using L for ln(n)
    print("The formula for P(n) is:")
    print(f"P(n) = ({coeff_L2_n2}*L^2 {coeff_L1_n2:+}*L + {coeff_L0_n2}) / ({denominator_n2}*n^2) + ({coeff_L3_n3}*L^3 {coeff_L2_n3:+}*L^2 + {coeff_L1_n3}*L) / ({denominator_n3}*n^3)")
    print("where L = ln(n).")

generate_formula()