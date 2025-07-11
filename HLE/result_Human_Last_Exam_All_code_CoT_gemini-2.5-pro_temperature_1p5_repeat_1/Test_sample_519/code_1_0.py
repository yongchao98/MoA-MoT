def solve_cfg_properties():
    """
    This function calculates the properties of the three given
    categories fibered in groupoids and prints them in the specified format.
    """

    # --- Properties for X1 = Hilb^11(A^3) ---
    # Type: S (Scheme)
    # Separated: s
    # Not universally closed, Not irreducible.
    # Dimension: d * n
    d1_degree = 11
    d1_affine_space_dim = 3
    dim1 = d1_degree * d1_affine_space_dim
    x1_profile = f"[S,s,{dim1}]"

    # --- Properties for X2 = [(A^4 \ V(xy-zw))/C*] ---
    # Type: S (Scheme, since action is free)
    # Separated: s
    # Irreducible: irr
    # Not universally closed.
    # Dimension: dim(space) - dim(group)
    dim_total_space = 4
    dim_group = 1
    dim2 = dim_total_space - dim_group
    x2_profile = f"[S,s,irr,{dim2}]"

    # --- Properties for X3 = Pic(C_0) for g=7 curve C_0 ---
    # Type: A (Algebraic Stack)
    # Separated: s
    # Not universally closed, Not irreducible.
    # Dimension: g - 1
    genus = 7
    dim3 = genus - 1
    x3_profile = f"[A,s,{dim3}]"

    # --- Final formatted output ---
    final_answer_string = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer_string)

solve_cfg_properties()