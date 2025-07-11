def solve_dessin_questions():
    """
    This function prints the solution to the multipart question about dessins d'enfants.
    """
    # Based on the step-by-step analysis:
    # (a) The quasiprimitive action on faces forces the socle N to be transitive,
    #     making the quotient dessin D_N unicellular. So, Yes.
    # (b) Smooth coverings are possible for types HA and AS, but not for type TW.
    # (c) Since smooth coverings are never possible for type TW, the statement for l <= 5 is true.
    
    answer_a = "Yes"
    answer_b = "HA, AS"
    answer_c = "True"
    
    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_dessin_questions()