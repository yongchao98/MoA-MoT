def solve():
    """
    This function determines and prints the Linear Logic formula for F(r,z,nz).
    """

    # The formula F(r,z,nz) must encode the choice in the JZDEC instruction.
    # The choice is represented by the additive conjunction '&'.

    # Branch 1: The counter is non-zero (x_r > 0).
    # Action: Decrement counter r and jump to state nz.
    # Encoding: Consume one C_r resource and produce S_nz.
    # This is represented by the linear implication C_r --o S_nz.
    # where --o is the multiset rewrite operator (multilinear implication).
    nonzero_branch = "C_r \\multimap S_{nz}"

    # Branch 2: The counter is zero (x_r = 0).
    # Action: Jump to state z.
    # Encoding: This branch must only succeed if no C_r is present.
    # We use the provided Z_r literal for the zero-test.
    # The formula introduces S_z for the new state and Z_r to perform the test.
    # The axioms in Delta ensure that Z_r will cause the proof to fail
    # if any C_r is present, and to succeed otherwise (by being consumed via Z_r --o 1).
    zero_branch = "S_z \\otimes Z_r"

    # The final formula combines these two mutually exclusive branches with '&'.
    # A correct derivation will be forced to choose the branch that matches
    # the state of the counter resources.
    F_formula = f"({zero_branch}) \\& ({nonzero_branch})"

    print("The formula F(r,z,nz) that correctly encodes the JZDEC instruction is:")
    print(F_formula)

solve()