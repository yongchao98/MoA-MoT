def solve_dessin_question():
    """
    Solves the theoretical question about dessin d'enfants and their automorphism groups.
    The reasoning is as follows:
    (a) Yes. G's quasiprimitive action on faces implies its socle N is transitive on faces.
        The orbits of N on the faces are the faces of the quotient D_N.
        Transitivity means there is one orbit, so D_N has one face (is unicellular).
    (b) HA, AS. For a group G to be the automorphism group of a regular dessin, it must be
        a (2, m, k)-group. Groups of type HA (affine) and AS (almost simple) are known to occur.
        Groups of type TW (wreath product T wr Z_l with T non-abelian simple, l>1) are known
        *not* to be (2, m, k)-groups, by a theorem of Conder.
    (c) True. As argued for (b), a group of type TW (l>1) cannot be the automorphism group of a regular
        dessin. The premise of the if-statement is thus false under the problem's conditions.
        An implication with a false premise is always true.
    """
    
    # Combine the answers into the specified format
    answer_a = "Yes"
    answer_b = "HA, AS"
    answer_c = "True"
    
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    
    print(final_answer)

solve_dessin_question()