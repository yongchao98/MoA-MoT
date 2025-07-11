def solve_group_theory_questions():
    """
    This function provides the answers to the theoretical questions about
    free groups and formal languages.
    """
    
    # (a) Is alpha(K) context-free for every rational subset K?
    # Yes, this is a known result in geometric group theory.
    answer_a = "Yes"
    
    # (b) If alpha(K) is context-free, does this imply that Geo(alpha(K)) is also context-free?
    # Yes, this follows from standard definitions and the closure of CFLs
    # under intersection with regular languages.
    answer_b = "Yes"
    
    # (c) If K is context-free, must Geo(alpha(K)) be context-free?
    # No, the conjugacy closure of a context-free language is not necessarily
    # context-free. A counterexample can be constructed.
    answer_c = "No"

    print(f"({answer_a}); ({answer_b}); ({answer_c})")

solve_group_theory_questions()