import sympy

def compute_homology():
    """
    This function computes the structure of the abelian group Z_2 x Z_3
    by finding the Smith Normal Form of its relation matrix.
    This corresponds to the first homology group H_1(X, Z).
    """
    # The abelian group G = Z_2 x Z_3 has a presentation <a, b | 2a=0, 3b=0>.
    # The relation matrix for this presentation is:
    relation_matrix = sympy.Matrix([[2, 0], [0, 3]])

    # The structure of the abelian group is given by the invariant factors,
    # which are the diagonal entries of the Smith Normal Form of the relation matrix.
    snf_matrix = relation_matrix.smith_normal_form()
    invariant_factors = [d for d in snf_matrix.diagonal() if d != 1]

    # The group is the direct sum of cyclic groups of these orders.
    # Z_1 is the trivial group, so we ignore factors of 1.
    group_structure = " x ".join([f"Z_{d}" for d in invariant_factors])

    # By the Chinese Remainder Theorem, Z_2 x Z_3 is isomorphic to Z_6.
    # The final equation is Z_2 x Z_3 = Z_6. The numbers involved are 2, 3, and 6.
    print("The final isomorphism is Z_2 x Z_3 = Z_6.")
    print("The numbers in this equation are:")
    print(2)
    print(3)
    final_order = 1
    for factor in invariant_factors:
        final_order *= factor
    print(final_order)

compute_homology()