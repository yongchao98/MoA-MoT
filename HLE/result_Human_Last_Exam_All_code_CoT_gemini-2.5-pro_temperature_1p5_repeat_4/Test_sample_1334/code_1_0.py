def solve_dessin_theory_questions():
    """
    This function provides the answers to the three theoretical questions
    about dessins d'enfants and their associated group theory.
    """

    # For part (a), the definition of a smooth covering implies that the covering
    # group N must act without fixed points on the vertices and faces, which is
    # the definition of a semi-regular action.
    answer_a = "True"

    # For part (b), the definition of a quasi-simple group is that it is a
    # perfect group (G = G') such that G/Z(G) is simple. The statement
    # is therefore true by definition.
    answer_b = "True"

    # For part (c), the conditions together imply that G is a quasiprimitive
    # group with a minimal normal subgroup N that acts regularly on the faces of the dessin.
    # The classification of such groups gives two possibilities: Holomorph Affine (HA) if N
    # is abelian, and Twisted Wreath (TW) if N is non-abelian.
    # However, a known theorem states that a group of type TW cannot be the automorphism
    # group of a regular dessin. This eliminates the TW case, leaving only HA.
    answer_c = "HA"

    # The final answer is formatted as requested.
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

solve_dessin_theory_questions()