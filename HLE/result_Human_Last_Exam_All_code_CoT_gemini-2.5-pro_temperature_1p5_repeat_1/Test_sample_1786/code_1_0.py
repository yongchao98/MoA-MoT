def generate_jzdec_formula():
    """
    This function constructs and prints the Linear Logic formula F(r, z, nz)
    for a JZDEC instruction in a Minsky machine.
    """
    
    # We will use symbolic variables r, z, and nz as per the problem description.
    register_var = 'r'
    zero_label_var = 'z'
    nonzero_label_var = 'nz'

    # The formula has two parts, corresponding to the two cases of JZDEC,
    # combined with the '&' (With) connective for external choice.

    # 1. Zero Branch: If counter 'r' is 0.
    # The state transitions to 'z'. In our logic, S_l is replaced by S_z.
    zero_branch = f"S_{{{zero_label_var}}}"

    # 2. Non-Zero Branch: If counter 'r' is > 0.
    # The state transitions to 'nz' and the counter 'r' is decremented.
    # This is encoded as an implication that consumes a C_r resource and produces S_nz.
    # The symbol '⊸' represents linear implication.
    nonzero_branch = f"(C_{{{register_var}}} ⊸ S_{{{nonzero_label_var}}})"
    
    # The final formula combines both branches. The proof context (presence or absence
    # of a C_r resource) will determine which branch can be successfully followed.
    final_formula = f"{zero_branch} & {nonzero_branch}"

    print("The formula for F(r, z, nz) that models the JZDEC instruction is:")
    print(final_formula)
    
    # To satisfy the prompt's instruction to "output each number",
    # here is a concrete example with r=2, z=4, nz=5.
    r, z, nz = 2, 4, 5
    example_formula = f"S_{{{z}}} & (C_{{{r}}} ⊸ S_{{{nz}}})"
    print(f"\nFor example, if r={r}, z={z}, and nz={nz}, the formula F({r},{z},{nz}) would be:")
    print(example_formula)


# Execute the function to print the formula.
generate_jzdec_formula()