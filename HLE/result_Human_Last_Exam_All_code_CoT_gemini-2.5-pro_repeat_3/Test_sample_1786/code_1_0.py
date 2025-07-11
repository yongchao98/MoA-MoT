def solve_linear_logic_formula():
    """
    This function determines and prints the formula for F(r, z, nz).
    """
    # The formula F(r, z, nz) is a conjunction of two cases:
    # 1. The Zero Case: Produces the new state literal S_z and a zero-test token Z_r.
    #    Represented by: (S_z ⊗ Z_r)
    zero_branch = "(S_z ⊗ Z_r)"

    # 2. The Non-Zero Case: Consumes one unit of counter r (C_r) and produces the new state literal S_nz.
    #    Represented by: (C_r ⊸ S_nz)
    #    Note: Using '⊸' for linear implication 'multimap'.
    nonzero_branch = "(C_r ⊸ S_nz)"

    # The two branches are combined with the additive conjunction '&' (with).
    # This allows the proof to choose a branch, but only the correct branch
    # (based on the context of the counters) will lead to a successful derivation.
    final_formula = f"{zero_branch} & {nonzero_branch}"

    # The variables r, z, and nz from the problem description are symbolic.
    # The prompt asks to output numbers in the final equation.
    # Since r, z, and nz are symbols, not numbers, we will print the formula
    # as a string, highlighting the symbolic components as requested.
    # We will "output" them by ensuring they are clearly part of the printed formula.
    # For example, if r=1, z=5, nz=6, the formula would be:
    # F(1, 5, 6) = (S_5 ⊗ Z_1) & (C_1 ⊸ S_6)
    # We will print the general symbolic formula.

    print("The formula F(r, z, nz) that encodes the JZDEC instruction is:")
    print(f"F(r, z, nz) = {final_formula}")
    print("\nWhere:")
    print("  - 'r' is the index of the counter being tested.")
    print("  - 'z' is the state label for the zero branch.")
    print("  - 'nz' is the state label for the non-zero branch.")
    print("  - '⊗' is the multiplicative conjunction (tensor).")
    print("  - '&' is the additive conjunction (with).")
    print("  - '⊸' is linear implication (lollipop).")

solve_linear_logic_formula()