def solve_cfg_properties():
    """
    This function analyzes three categories fibered in groupoids and determines their properties.
    """
    # Properties for X1
    x1_type = "S"
    x1_props = ["s"]
    x1_dim = 33
    x1_profile = f"[{','.join([x1_type] + x1_props + [str(x1_dim)])}]"

    # Properties for X2
    x2_type = "DM"
    x2_props = ["s", "irr"]
    x2_dim = 3
    x2_profile = f"[{','.join([x2_type] + x2_props + [str(x2_dim)])}]"

    # Properties for X3
    x3_type = "A"
    x3_props = ["s"]
    x3_dim = 6
    x3_profile = f"[{','.join([x3_type] + x3_props + [str(x3_dim)])}]"

    # Final combined answer
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer)

solve_cfg_properties()