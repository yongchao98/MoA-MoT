import collections

def print_formula(r, z, nz):
    """
    Constructs and prints the Linear Logic formula F(r, z, nz).

    Args:
      r: The register index (integer).
      z: The state label for the zero case (string).
      nz: The state label for the non-zero case (string).
    """
    
    # Using unicode for logic symbols for better readability
    tensor = "\u2297"  # ⊗
    multimap = "\u22B8" # ⊸
    with_op = "&"      # &

    # Non-zero case: Decrement counter r and transition to state nz.
    # This requires consuming a C_r resource and producing an S_nz resource.
    # Formula: C_r ⊸ S_nz
    non_zero_branch = f"C_{r} {multimap} S_{nz}"

    # Zero case: Test that counter r is zero and transition to state z.
    # This produces the new state S_z and a zero-test token Z_r.
    # This branch will fail during proof search if a C_r token is present.
    # Formula: S_z ⊗ Z_r
    zero_branch = f"S_{z} {tensor} Z_{r}"

    # F(r,z,nz) is the additive conjunction ("with") of the two branches.
    # This provides an internal choice that is determined by the context (presence of C_r).
    final_formula = f"({non_zero_branch}) {with_op} ({zero_branch})"

    print("The formula F(r,z,nz) is:")
    print(final_formula)
    
    # Example usage:
    # Let's say we have an instruction JZDEC(1, "state_z", "state_nz")
    # r=1, z="state_z", nz="state_nz"
    print("\nFor example, F(1, z, nz) would be:")
    
    example_r = 1
    example_z = 'z'
    example_nz = 'nz'
    
    example_non_zero = f"C_{example_r} {multimap} S_{example_nz}"
    example_zero = f"S_{example_z} {tensor} Z_{example_r}"
    example_formula = f"({example_non_zero}) {with_op} ({example_zero})"
    print(example_formula)


# You can change these values to see the formula for a different instruction.
register_index = 1
zero_state_label = "z"
non_zero_state_label = "nz"

print_formula(register_index, zero_state_label, non_zero_state_label)
