def generate_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for the JZDEC instruction.

    Args:
        r (int or str): The register index (1-based).
        z (str): The state label for the zero case.
        nz (str): The state label for the non-zero case.
    """
    
    # Formula for the case where counter r is 0.
    # It introduces the test literal Z_r and the new state literal S_z.
    # S_z tensor Z_r
    zero_case_formula = f"S_{z} \u2297 Z_{r}"

    # Formula for the case where counter r > 0.
    # It consumes a C_r literal and produces the new state literal S_nz.
    # C_r multimap S_nz
    nonzero_case_formula = f"C_{r} \u27f8 S_{nz}"

    # The full formula uses the additive conjunction '&' to represent the
    # conditional choice. The prover must pick the branch that is provable
    # with the current resources.
    # (S_z tensor Z_r) with (C_r multimap S_nz)
    F_r_z_nz = f"({zero_case_formula}) & ({nonzero_case_formula})"

    print("The formula F(r, z, nz) is:")
    print(F_r_z_nz)
    print("\nExplanation of the components:")
    print(f"1. Non-zero branch '{nonzero_case_formula}': If counter {r} is positive, it has at least one C_{r} resource.")
    print(f"   This formula consumes one C_{r} (decrementing the counter) and produces the new state S_{nz}.")
    print(f"2. Zero branch '{zero_case_formula}': If counter {r} is zero, there is no C_{r} resource.")
    print(f"   This formula introduces S_{z} for the new state and Z_{r} to begin the zero-test.")
    print(f"   The axioms in \u0394 will use Z_{r} to ensure no C_{r} exists. This branch fails if C_{r} is present.")

# Example usage with r=1, z='q_z', nz='q_nz'
generate_formula(r=1, z='z', nz='nz')