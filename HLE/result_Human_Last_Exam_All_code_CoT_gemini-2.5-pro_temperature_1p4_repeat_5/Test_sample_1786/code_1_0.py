def solve_linear_logic_formula():
    """
    This script determines and prints the Linear Logic formula for F(r, z, nz)
    that correctly models a Minsky machine's JZDEC instruction.
    """

    # The variables r, z, and nz are symbolic.
    # r represents the register index (a number).
    # z represents the state label for the zero branch.
    # nz represents the state label for the non-zero branch.
    r = "r"
    z = "z"
    nz = "nz"

    # The formula F(r, z, nz) is a choice between two cases, modeled by the '&' (with) connective.

    # Case 1: The counter x_r is zero. The new state is 'z'.
    # The formula for this branch simply introduces the new state literal.
    zero_branch_formula = f"S_{z}"

    # Case 2: The counter x_r is non-zero. The state becomes 'nz' and counter r is decremented.
    # This corresponds to consuming a C_r resource and producing an S_nz resource.
    # This is encoded by a linear implication, where '-o' represents the 'multimap' symbol (⊸).
    nonzero_branch_formula = f"(C_{r} -o S_{nz})"

    # The final formula combines both cases with the '&' connective.
    final_formula = f"{zero_branch_formula} & {nonzero_branch_formula}"

    print("The appropriate formula for F(r, z, nz) is:")
    print(final_formula)
    
    print("\nIn equation form, this is:")
    # This prints the final "equation" for F, showing the symbolic numbers/labels r, z, and nz.
    print(f"F({r}, {z}, {nz}) = {final_formula}")

solve_linear_logic_formula()