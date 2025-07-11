def solve_cfgs():
    """
    This function calculates the properties of the three given CFGs and prints the result.
    """

    # Properties for X_1: Hilb^11(A^3)
    # Type: Scheme (S)
    # Separated: yes (s)
    # Universally Closed: no
    # Irreducible: no
    # Dimension: n*d for the main component of Hilb^d(A^n)
    d1 = 11  # degree
    n1 = 3   # dimension of affine space
    dim1 = n1 * d1
    x1_profile = f"[S, s, {dim1}]"

    # Properties for X_2: [(A^4 \ V(xy-zw)) / C*]
    # Type: Deligne-Mumford stack (DM), since stabilizers are finite
    # Separated: yes (s)
    # Universally Closed: no
    # Irreducible: yes (irr)
    # Dimension: dim(space) - dim(group)
    dim_space2 = 4
    dim_group2 = 1
    dim2 = dim_space2 - dim_group2
    x2_profile = f"[DM, s, irr, {dim2}]"

    # Properties for X_3: Pic_{C_0} for a genus 7 curve C_0
    # Type: Algebraic stack (A), since stabilizers are C*
    # Separated: yes (s)
    # Universally Closed: no
    # Irreducible: no
    # Dimension: dim(Jacobian) + dim(stabilizer) = g + 1
    g3 = 7  # genus
    dim_stab3 = 1 # dim(C*)
    dim3 = g3 + dim_stab3
    x3_profile = f"[A, s, {dim3}]"

    # Combine and print the results
    final_answer = f"{x1_profile} {x2_profile} {x3_profile}"
    print(final_answer)

solve_cfgs()