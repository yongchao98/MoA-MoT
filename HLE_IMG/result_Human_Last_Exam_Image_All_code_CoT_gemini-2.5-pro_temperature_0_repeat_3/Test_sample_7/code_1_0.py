def solve_formula():
    """
    This function determines and prints the simple formula for the second voltage
    plateau of a graphite anode during lithium intercalation.
    """

    # The second voltage plateau corresponds to the formation of Stage 2 material
    # from Stage 3 material. In the notation V_k for the plateau forming stage k,
    # the second plateau is V_2.
    k = 2

    # The stage being consumed is the next higher stage, k+1.
    k_plus_1 = k + 1

    # A simple model for the voltage of the plateau between stage k+1 and stage k
    # is given by the difference in their respective chemical potentials (μ_k and μ_{k+1}),
    # divided by the elementary charge (e).
    # The formula is V_k = (μ_k - μ_{k+1}) / e.

    # We construct the final formula string for k=2.
    # The numbers in the equation are the indices k and k+1.
    final_formula = f"(μ_{k} - μ_{k_plus_1}) / e"

    print("A simple formula that approximates the second plateau is:")
    print(final_formula)

solve_formula()