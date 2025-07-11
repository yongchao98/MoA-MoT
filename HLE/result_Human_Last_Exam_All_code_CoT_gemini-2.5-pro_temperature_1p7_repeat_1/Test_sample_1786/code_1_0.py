def print_formula(r, z, nz):
    """
    Prints the Linear Logic formula F(r, z, nz).

    Args:
      r: The index of the register (1-based).
      z: The instruction label for the zero case.
      nz: The instruction label for the non-zero case.
    """
    # The formula F(r,z,nz) is S_z & (C_r -o S_nz)
    # where '&' is the additive conjunction (With)
    # and '-o' is the linear implication (lollipop)

    # Let's output each part of the formula clearly.
    # The state for the zero branch:
    case_zero = f"S_{z}"

    # The formula for the non-zero, decrementing branch:
    case_nonzero_antecedent = f"C_{r}"
    case_nonzero_consequent = f"S_{nz}"
    case_nonzero = f"({case_nonzero_antecedent} \u22B8 {case_nonzero_consequent})" # Using ⊸ for linear implication

    # The full formula F combines these with the '&' connective:
    # We will use the ampersand symbol '&' to represent the 'With' connective.
    final_formula = f"{case_zero} & {case_nonzero}"

    print(f"The formula F(r, z, nz) for the JZDEC instruction is:")
    print(final_formula)
    print("\nBreaking it down:")
    print(f"  - Instruction label for 'zero' case: {z}")
    print(f"  - Instruction label for 'non-zero' case: {nz}")
    print(f"  - Register index: {r}")
    print("\n  - Component for zero case (state becomes z, counters unchanged):")
    print(f"    S_{z}")
    print("\n  - Component for non-zero case (decrement r, state becomes nz):")
    print(f"    C_{r} \u22B8 S_{nz}")
    print("\nThese components are combined with the 'With' connective (&), representing a choice:")
    print(f"F({r}, {z}, {nz}) = S_{z} & (C_{r} \u22B8 S_{nz})")


# Example usage with some symbolic values
# Let's take r=1, z='z_label', nz='nz_label'
print_formula(r=1, z='z', nz='nz')