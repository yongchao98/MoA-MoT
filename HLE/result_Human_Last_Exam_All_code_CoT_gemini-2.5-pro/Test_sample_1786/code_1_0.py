def generate_jzdec_formula(r, z, nz):
    """
    This function generates the Linear Logic formula F(r, z, nz)
    for the JZDEC instruction.

    Args:
        r (str or int): The register to test.
        z (str): The label for the zero case.
        nz (str): The label for the non-zero case.
    """
    # The formula combines the two cases using the '&' (With) connective.
    # Case 1 (Zero): S_z
    #   If the counter is 0, the machine transitions to state z.
    # Case 2 (Non-Zero): (C_r -o S_nz)
    #   If the counter is non-zero, it consumes a C_r resource (decrement)
    #   and transitions to state nz. The '-o' is linear implication.
    
    # The & connective on the left of a sequent allows the prover to choose
    # which branch to follow. If the counter r is 0, the (C_r -o S_nz) branch
    # is impossible, forcing the correct choice. If the counter is > 0,
    # choosing the S_z branch corresponds to an incorrect machine transition,
    # which will not lead to a valid proof of acceptance.

    formula = f"S_{z} & (C_{r} -o S_{nz})"
    print("The formula F(r, z, nz) for the JZDEC instruction is:")
    print(formula)

# Example usage with symbolic parameters r, z, nz as in the problem description.
generate_jzdec_formula('r', 'z', 'nz')

# The final answer as a raw string is S_z & (C_r -o S_{nz})
# The problem asks to output the answer in a specific format at the end.
# The format seems to be expecting the formula itself.
final_answer = "S_z & (C_r -o S_{nz})"
print(f"\n<<<__{final_answer}__>>>")