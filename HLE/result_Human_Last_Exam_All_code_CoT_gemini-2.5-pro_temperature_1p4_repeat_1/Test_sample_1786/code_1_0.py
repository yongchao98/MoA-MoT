def generate_and_print_formula():
    """
    This function constructs and prints the Linear Logic formula for the JZDEC instruction.
    It first prints the symbolic formula and then an example with concrete numbers.
    """

    # Define symbolic components of the formula
    z_r = "Z_r"
    s_z = "S_z"
    c_r = "C_r"
    s_nz = "S_nz"

    # Use unicode characters for Linear Logic connectives for readability
    tensor = "\u2297"  # ⊗ (otimes)
    lolli = "\u22B8"   # ⊸ (multimap)
    with_op = "&"      # & (with)

    # Construct the two branches of the internal choice
    zero_branch = f"({z_r} {lolli} {s_z})"
    nonzero_branch = f"({c_r} {tensor} {z_r} {lolli} {s_nz})"

    # Assemble the final formula F(r, z, nz)
    symbolic_formula = f"{z_r} {tensor} ({zero_branch} {with_op} {nonzero_branch})"

    print("The symbolic formula for F(r, z, nz) is:")
    print(symbolic_formula)
    print("-" * 20)

    # Generate an example as requested, with numbers for r, z, and nz
    r_ex, z_ex, nz_ex = 1, 5, 6
    z_r_ex = f"Z_{r_ex}"
    s_z_ex = f"S_{z_ex}"
    c_r_ex = f"C_{r_ex}"
    s_nz_ex = f"S_{nz_ex}"

    zero_branch_ex = f"({z_r_ex} {lolli} {s_z_ex})"
    nonzero_branch_ex = f"({c_r_ex} {tensor} {z_r_ex} {lolli} {s_nz_ex})"
    
    example_formula = f"{z_r_ex} {tensor} ({zero_branch_ex} {with_op} {nonzero_branch_ex})"

    print(f"An example formula for F(r={r_ex}, z={z_ex}, nz={nz_ex}), showing each number:")
    print(example_formula)
    
    # The final answer wrapped as requested
    final_answer = f"Z_r \u2297 ((Z_r \u22B8 S_z) & (C_r \u2297 Z_r \u22B8 S_{nz}))"
    print("\n<<<" + final_answer + ">>>")

# Execute the function to perform the task
generate_and_print_formula()