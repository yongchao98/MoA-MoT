def generate_jzdec_formula(r, z, nz):
    """
    Generates and prints the Linear Logic formula for a JZDEC instruction.

    Args:
        r (int): The register to test (e.g., 1, 2, ...).
        z (str or int): The label for the zero branch.
        nz (str or int): The label for the non-zero branch.
    """
    
    # Using common text representations for Linear Logic connectives:
    # ⊗ (tensor): a multiplicative conjunction
    # ⊸ (multimap/lollipop): linear implication
    # & (with): an additive conjunction
    
    # Construct the two main branches of the formula
    zero_branch = f"(S_{z} ⊗ Z_{r})"
    nonzero_branch = f"(C_{r} ⊸ S_{nz})"
    
    # Combine them with the 'with' connective
    final_formula = f"{zero_branch} & {nonzero_branch}"
    
    print("The appropriate formula F(r, z, nz) is:")
    print(final_formula)
    
    # To comply with the instruction "output each number in the final equation",
    # we show the symbolic representation again with placeholders.
    print("\nSymbolic representation:")
    print(f"F({r}, {z}, {nz}) = (S_{z} ⊗ Z_{r}) & (C_{r} ⊸ S_{nz})")

if __name__ == '__main__':
    # Example Usage:
    # Consider the instruction for register r=1,
    # with a jump to label L3 if zero, and to L4 if non-zero.
    example_r = 1
    example_z = "L3"
    example_nz = "L4"
    
    print(f"Example for instruction JZDEC(r={example_r}, z='{example_z}', nz='{example_nz}'):")
    print("-" * 50)
    generate_jzdec_formula(example_r, example_z, example_nz)
    print("-" * 50)
    
    # Example with numeric labels as requested
    print("\nExample with numeric labels JZDEC(r=2, z=5, nz=6):")
    print("-" * 50)
    generate_jzdec_formula(2, 5, 6)
    print("-" * 50)
    print("<<< (S_z ⊗ Z_r) & (C_r ⊸ S_{nz}) >>>")