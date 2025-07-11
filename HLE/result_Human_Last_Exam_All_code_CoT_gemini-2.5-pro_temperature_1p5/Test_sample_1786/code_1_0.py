def generate_jzdec_formula(r, z, nz):
    """
    Generates the Linear Logic formula F(r, z, nz) for a JZDEC instruction.

    Args:
        r: The register index (string or int).
        z: The state label for the zero case (string or int).
        nz: The state label for the non-zero case (string or int).

    Returns:
        A string representing the formula.
    """
    # Using Unicode for logical operators for better readability
    tensor = "\u2297"  # ⊗
    multimap = "\u22b8" # ⊸
    with_op = "&"      # &

    # Build the two branches of the formula
    # Branch 1: The case where counter r is zero.
    # It produces the new state S_z and a 'zero-checker' atom Z_r.
    zero_branch = f"(S_{z} {tensor} Z_{r})"

    # Branch 2: The case where counter r is non-zero.
    # It consumes a C_r and produces the new state S_nz.
    nonzero_branch = f"(C_{r} {multimap} S_{nz})"

    # Combine the branches with the additive conjunction '&'
    formula = f"{zero_branch} {with_op} {nonzero_branch}"
    
    return formula

def main():
    """
    Main function to define symbolic variables and print the final formula.
    """
    # Symbolic placeholders for the instruction parameters
    r_sym = 'r'
    z_sym = 'z'
    nz_sym = 'nz'

    # Generate the formula
    final_formula = generate_jzdec_formula(r_sym, z_sym, nz_sym)
    
    # Print the final result in a clear format
    print("The appropriate formula F(r, z, nz) is:")
    print(f"F({r_sym}, {z_sym}, {nz_sym}) = {final_formula}")
    
    print("\nWhere the components are:")
    print(f"  - S_{z_sym}: Literal for the 'zero' destination state '{z_sym}'.")
    print(f"  - S_{nz_sym}: Literal for the 'non-zero' destination state '{nz_sym}'.")
    print(f"  - C_{r_sym}: Literal for a unit in counter '{r_sym}'.")
    print(f"  - Z_{r_sym}: Literal used to test if counter '{r_sym}' is zero.")
    print(f"  - {tensor} (tensor): Multiplicative conjunction.")
    print(f"  - & (with): Additive conjunction, representing an external choice.")
    print(f"  - {multimap} (multimap): Linear implication.")

if __name__ == "__main__":
    main()
