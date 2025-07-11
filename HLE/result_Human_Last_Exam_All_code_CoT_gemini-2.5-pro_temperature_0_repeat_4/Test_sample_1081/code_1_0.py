def get_number_of_tori():
    """
    This function calculates and prints the formula for the exact number of
    F_q-rational maximal tori of a reductive group of type E_8.
    """
    # Step 1: Define the parameters for the E_8 group type.
    # The number of positive roots for E_8.
    num_positive_roots = 120
    # The rank of E_8.
    rank = 8
    # The degrees of the fundamental invariants for the Weyl group of type E_8.
    degrees = [2, 8, 12, 14, 18, 20, 24, 30]

    # Step 2: Explain the theoretical basis for the calculation.
    print("A reductive group G of type E_8 over a finite field F_q is always split.")
    print("The number of F_q-rational maximal tori is given by the ratio |G(F_q)| / |T_0(F_q)|,")
    print("where T_0 is a split maximal torus.")
    print("")

    # Step 3: Present the formula for the order of the group G(F_q).
    order_G_terms = [f"q^{num_positive_roots}"]
    for d in degrees:
        order_G_terms.append(f"(q^{d} - 1)")
    print(f"The order of the group is |G(F_q)| = {' * '.join(order_G_terms)}.")

    # Step 4: Present the formula for the order of the split torus T_0(F_q).
    print(f"The order of the split torus is |T_0(F_q)| = (q - 1)^{rank}.")
    print("")

    # Step 5: Construct and print the final formula for the number of tori.
    print("The exact number of F_q-rational maximal tori is the ratio, which simplifies to:")
    
    final_formula_terms = [f"q^{num_positive_roots}"]
    for d in degrees:
        # The term (q^d - 1) / (q - 1) is the sum of a geometric series.
        # We present it in a way that makes the degree 'd' explicit.
        final_formula_terms.append(f"((q^{d} - 1)/(q - 1))")
    
    final_formula = " * ".join(final_formula_terms)
    print(final_formula)

get_number_of_tori()