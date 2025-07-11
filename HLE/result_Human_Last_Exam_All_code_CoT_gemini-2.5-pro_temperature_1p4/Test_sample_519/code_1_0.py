def solve_cfg_properties():
    """
    This function determines the properties of the three given CFGs
    and prints the formatted result.
    """

    # --- X1: Hilb_11(A^3) ---
    # Type: Scheme (S)
    # Properties: Not separated, not universally closed, not irreducible.
    # Dimension calculation: n * d = 3 * 11 = 33
    x1_type = "S"
    x1_dim = 33
    x1_properties = [x1_type, str(x1_dim)]
    x1_profile = f"[{','.join(x1_properties)}]"

    # --- X2: [ (A^4 \ V(xy-zw)) / C* ] ---
    # Type: Deligne-Mumford Stack (DM)
    # Properties: separated (s), irreducible (irr)
    # Dimension calculation: dim(A^4) - dim(C*) = 4 - 1 = 3
    x2_type = "DM"
    x2_properties = ["s", "irr"]
    x2_dim = 3
    x2_all_props = [x2_type] + x2_properties + [str(x2_dim)]
    x2_profile = f"[{','.join(x2_all_props)}]"

    # --- X3: Pic(C_0) for a genus 7 curve C_0 ---
    # Type: Scheme (S)
    # Properties: separated (s)
    # Dimension calculation: genus g = 7
    x3_type = "S"
    x3_properties = ["s"]
    x3_dim = 7
    x3_all_props = [x3_type] + x3_properties + [str(x3_dim)]
    x3_profile = f"[{','.join(x3_all_props)}]"

    # Combine the profiles into a single string
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"

    print(final_answer)

solve_cfg_properties()