def generate_jzdec_formula():
    """
    This function constructs the Linear Logic formula F(r, z, nz)
    for the JZDEC instruction and prints it.
    """
    # Symbolic placeholders for the instruction parameters
    r = 'r'   # The register index
    z = 'z'   # The state label for the zero branch
    nz = 'nz' # The state label for the non-zero branch

    # The formula for the zero-case branch.
    # It produces a "zero-test" token Z_r along with the new state S_z.
    # The helper axioms are set up so that Z_r can only be resolved
    # (i.e., removed from the sequent) if no C_r token is present.
    zero_branch_formula = f"(Z_{r} \\otimes S_{z})"

    # The formula for the non-zero-case branch.
    # It consumes one counter resource C_r and produces the new state S_nz.
    # This is a direct encoding of "decrement and transition" using linear implication.
    nonzero_branch_formula = f"(C_{r} \\multimap S_{{nz}})"

    # The complete formula for F(r, z, nz) combines both cases using the
    # additive conjunction '&' (with). This represents an external choice,
    # where the context (the presence or absence of C_r) determines
    # which of the two branches is provable.
    final_formula = f"{zero_branch_formula} & {nonzero_branch_formula}"

    print("The appropriate formula for F(r, z, nz) is:")
    print(final_formula)
    print("\nWhere the components of the formula are:")
    print(f"1. Zero Case ('z'): {zero_branch_formula}")
    print(f"2. Non-Zero Case ('nz'): {nonzero_branch_formula}")

generate_jzdec_formula()