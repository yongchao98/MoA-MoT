def generate_jzdec_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for a JZDEC instruction.

    Args:
        r (int): The index of the counter register to test.
        z (int or str): The label of the state to jump to if the counter is zero.
        nz (int or str): The label of the state to jump to if the counter is non-zero.
    """
    # Unicode characters for logic symbols for better readability
    tensor = "\u2297"    # ⊗
    lollipop = "\u22B8" # ⊸
    with_op = "&"       # &

    # Construct the two branches of the formula
    zero_branch = f"(S_{z} {tensor} Z_{r})"
    nonzero_branch = f"(C_{r} {lollipop} S_{nz})"

    # Combine them using the 'with' connective
    formula = f"{zero_branch} {with_op} {nonzero_branch}"

    # Print the final result
    print(f"The formula F({r}, {z}, {nz}) for the JZDEC instruction is:")
    print(formula)
    print("\nwhere the full instruction is encoded as: ")
    print(f"!(S_l {lollipop} ({formula}))")


# Example usage with r=1, z='z', nz='nz'
# These are symbolic but the function can take numbers as well.
generate_jzdec_formula(r=1, z="z", nz="nz")

# Example with numeric labels as used in some Minsky machine definitions
generate_jzdec_formula(r=2, z=5, nz=6)
