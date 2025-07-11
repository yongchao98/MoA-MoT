def print_lower_bound_formula():
    """
    This function prints the derived symbolic formula for the lower bound on R_n^*.
    """

    # Symbolic names for the parameters and distributions in the problem.
    Phi = "Φ"
    delta = "δ"
    N = "N"
    P0_n = "P_0^n"
    Pj_n_sum = "Σ_{j=1 to N} P_j^n"

    # The constant numbers appearing in the formula.
    numerator_constant = 1
    denominator_constant_1 = 2
    denominator_constant_2 = 2

    # Assemble the formula as a string.
    # The bound is: (Φ(δ/2) / 2) * (1 - d_TV(P_0^n, (1/N) * Σ P_j^n))
    formula = (f"R_n* >= ({Phi}({delta}/{denominator_constant_1}) / {denominator_constant_2}) * "
               f"({numerator_constant} - d_TV({P0_n}, (1/{N}) * {Pj_n_sum}))")

    print("The tightest lower bound that can be proven with the given information is based on a reduction from estimation to hypothesis testing.")
    print("The final formula is:")
    print(formula)

print_lower_bound_formula()
