def solve_cfg_properties():
    """
    This function analyzes the properties of three categories fibered in groupoids
    and prints the results in the specified format.
    """

    # --- Analysis of X1 ---
    # X1 is the Hilbert scheme of 11 points in A^3.
    # Type: Scheme (S)
    # Separated: yes (s)
    # Universally closed: no
    # Irreducible: no
    # Dimension = n * d
    n_x1 = 3
    d_x1 = 11
    dim_x1 = n_x1 * d_x1
    profile_x1 = f"[S,s,{dim_x1}]"

    # --- Analysis of X2 ---
    # X2 is an open subscheme of the weighted projective space P(1,4,2,3).
    # Type: Scheme (S)
    # Separated: yes (s)
    # Universally closed: no
    # Irreducible: yes (irr)
    # Dimension = dim(A^4) - dim(C*)
    dim_space_x2 = 4
    dim_group_x2 = 1
    dim_x2 = dim_space_x2 - dim_group_x2
    profile_x2 = f"[S,s,irr,{dim_x2}]"

    # --- Analysis of X3 ---
    # X3 is the Picard scheme of a genus 7 curve.
    # Type: Scheme (S)
    # Separated: yes (s)
    # Universally closed: no
    # Irreducible: no
    # Dimension = genus
    g_x3 = 7
    dim_x3 = g_x3
    profile_x3 = f"[S,s,{dim_x3}]"

    # --- Final Output ---
    final_result = f"{profile_x1} {profile_x2} {profile_x3}"
    print(final_result)

solve_cfg_properties()