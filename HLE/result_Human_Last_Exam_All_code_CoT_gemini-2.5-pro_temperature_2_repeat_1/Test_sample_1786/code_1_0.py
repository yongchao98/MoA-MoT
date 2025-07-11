def generate_formula(r, z, nz):
    """
    Generates the Linear Logic formula for a JZDEC instruction.

    Args:
      r: The index of the register to check (an integer).
      z: The instruction label for the zero case (a string).
      nz: The instruction label for the non-zero case (a string).

    Returns:
      A string representing the formula F(r, z, nz).
    """
    # The non-zero branch: consumes a C_r and produces S_nz.
    # C_r -o S_nz
    f_nonzero = f"(C_{r} \multimap S_{nz})"

    # The zero branch: produces S_z and the zero-test gadget Z_r.
    # S_z (tensor) Z_r
    f_zero = f"(S_{z} \otimes Z_{r})"

    # Combine them with external choice '&'
    # The full formula F(r, z, nz) is the conjunction of the two cases.
    final_formula = f"{f_zero} & {f_nonzero}"

    print("The formula F(r,z,nz) encoding the instruction l: JZDEC(r, z, nz) is:")
    print(final_formula)
    print("\nExplanation:")
    print("This formula presents an external choice (&) to the proof system.")
    print("1. The first branch is chosen if the proof can succeed with it:")
    print(f"   {f_zero}")
    print(f"   This branch introduces a Z_{r} gadget. The provided axioms for Z_{r} make this branch succeed only if counter r is zero. S_{z} represents the next state.")
    print("2. The second branch is chosen otherwise:")
    print(f"   {f_nonzero}")
    print(f"   This is a linear implication. It succeeds only if a C_{r} literal (representing a non-zero counter r) is present to be consumed. It produces S_{nz} for the next state.")

# Example usage:
# l: JZDEC(r=1, z='state_A', nz='state_B')
generate_formula(r=1, z='A', nz='B')

# <<<(S_{z} \otimes Z_{r}) & (C_{r} \multimap S_{nz})>>>