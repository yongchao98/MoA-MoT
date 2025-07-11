def print_minsky_jzdec_formula():
    """
    This function prints the Linear Logic formula for the JZDEC instruction.
    """

    # Symbolic representations of the Minsky machine components
    state_z = "S_z"
    state_nz = "S_nz"
    counter_r = "C_r"
    zero_test_r = "Z_r"

    # Unicode characters for Linear Logic operators
    tensor_op = "\u2297"      # Multiplicative Conjunction (tensor)
    multimap_op = "\u27f8"     # Linear Implication (multimap)
    with_op = "&"            # Additive Conjunction (with)

    # Formula for the zero case: If x_r = 0, transition to state z.
    # The token Z_r is produced to interact with the axioms in Delta
    # to verify the absence of C_r tokens.
    zero_case_formula = f"({state_z} {tensor_op} {zero_test_r})"

    # Formula for the non-zero case: If x_r > 0, decrement x_r and transition to state nz.
    # This consumes one C_r token and produces a state token S_nz.
    non_zero_case_formula = f"({counter_r} {multimap_op} {state_nz})"

    # The two cases are combined with the '&' connective.
    # This creates a conditional branch where the context (the set of available resources)
    # determines which path can be successfully followed.
    final_formula = f"{zero_case_formula} {with_op} {non_zero_case_formula}"

    print("The appropriate formula for F(r,z,nz) is:")
    print(final_formula)

print_minsky_jzdec_formula()