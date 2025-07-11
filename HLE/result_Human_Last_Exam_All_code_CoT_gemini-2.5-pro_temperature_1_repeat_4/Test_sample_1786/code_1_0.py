def generate_jzdec_formula():
    """
    This function generates and prints the Linear Logic formula F(r, z, nz)
    for encoding the JZDEC instruction of a Minsky machine.
    """
    # Symbolic placeholders for the instruction parameters
    # r: the register index (a number, represented symbolically)
    # z: the instruction label for the zero case
    # nz: the instruction label for the non-zero case
    r = "r"
    z = "z"
    nz = "nz"

    # Unicode representations of Linear Logic operators
    tensor = "\u2297"    # ⊗ (otimes)
    multimap = "\u22B8"  # ⊸ (multimap)
    with_op = "&"        # & (with)

    # Formula for the case where counter r is zero.
    # It produces the new state S_z and a zero-test literal Z_r.
    # This branch is only provable if counter r is actually zero.
    zero_case = f"(S_{z} {tensor} Z_{r})"

    # Formula for the case where counter r is non-zero.
    # It consumes one unit of counter r (C_r) and produces the new state S_nz.
    # This branch is only provable if counter r is non-zero.
    non_zero_case = f"(C_{r} {multimap} S_{nz})"

    # The full formula F(r,z,nz) combines the two cases. The `&` operator
    # presents an external choice, and the context (the state of the counters)
    # determines which branch can be successfully proven.
    final_formula = f"{zero_case} {with_op} {non_zero_case}"

    print("The formula F(r, z, nz) that encodes the JZDEC instruction is:")
    print(final_formula)

if __name__ == '__main__':
    generate_jzdec_formula()