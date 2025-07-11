def generate_formula():
    """
    This function generates and prints the formula F(r, z, nz).
    """
    # Variables representing parts of the Minsky machine instruction
    r = "r"   # The register to test
    z = "z"   # The state to jump to if the register is zero
    nz = "nz" # The state to jump to if the register is non-zero

    # Literals used in the Linear Logic encoding
    S_z = f"S_{z}"
    S_nz = f"S_{nz}"
    C_r = f"C_{r}"
    Z_r = f"Z_{r}"

    # The two branches of the JZDEC instruction
    # Non-zero branch: consume a Cr resource and produce the new state S_nz
    nonzero_branch = f"({C_r} \multimap {S_nz})"
    
    # Zero branch: introduce the tester atom Z_r and the new state S_z
    # The Z_r atom will ensure this branch only succeeds if no Cr resources are present
    zero_branch = f"({Z_r} \otimes {S_z})"

    # The final formula F combines the two branches with the additive conjunction '&'
    # This represents the deterministic choice the machine makes based on its state
    F_formula = f"{zero_branch} & {nonzero_branch}"

    print("The formula F(r, z, nz) that encodes the JZDEC instruction is:")
    print(F_formula)

generate_formula()