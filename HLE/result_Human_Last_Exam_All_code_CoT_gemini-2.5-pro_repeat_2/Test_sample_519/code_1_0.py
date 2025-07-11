def solve_cfgs():
    """
    This function calculates the properties of the three given Categories Fibered in Groupoids (CFGs)
    and prints them in the specified format.
    """

    # --- X_1: Hilbert scheme of 11 points in A^3 ---
    # Type: Scheme (S)
    # Properties: separated (s)
    # Not universally closed, not irreducible.
    # Dimension calculation: n * d, for n=3 dimensions and d=11 points.
    dim_n_x1 = 3
    dim_d_x1 = 11
    dim_x1 = dim_n_x1 * dim_d_x1
    properties_x1 = ['S', 's', str(dim_x1)]
    profile_x1 = f"[{','.join(properties_x1)}]"

    # --- X_2: Quotient stack [(A^4 \ V(xy-zw))/C*] ---
    # Type: Deligne-Mumford stack (DM)
    # Properties: separated (s), irreducible (irr)
    # Not universally closed.
    # Dimension calculation: dim(ambient) - dim(group).
    dim_ambient_x2 = 4
    dim_group_x2 = 1
    dim_x2 = dim_ambient_x2 - dim_group_x2
    properties_x2 = ['DM', 's', 'irr', str(dim_x2)]
    profile_x2 = f"[{','.join(properties_x2)}]"

    # --- X_3: Picard stack of a genus 7 curve ---
    # Type: Algebraic stack (A)
    # Properties: separated (s)
    # Not universally closed, not irreducible.
    # Dimension calculation: genus + dim(stabilizer C*).
    genus_x3 = 7
    dim_stabilizer_x3 = 1
    dim_x3 = genus_x3 + dim_stabilizer_x3
    properties_x3 = ['A', 's', str(dim_x3)]
    profile_x3 = f"[{','.join(properties_x3)}]"

    # --- Final Answer ---
    # Combine the profiles and print the final result.
    final_answer = f"{profile_x1} {profile_x2} {profile_x3}"
    print(final_answer)

solve_cfgs()