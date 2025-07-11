def solve_complexity_questions():
    """
    Solves the theoretical complexity theory questions based on definitions.
    """

    # (a) If the problem of deciding whether S is not free is NP-hard, does that imply
    # the problem of deciding whether S is free is also NP-hard?
    # No. If "not free" is NP-hard, then "free" is co-NP-hard. An NP-hard problem
    # is not necessarily co-NP-hard, and vice-versa, unless NP = co-NP.
    answer_a = "No"

    # (b) If the problem of deciding whether S is not free is NP-complete, does that imply
    # the problem of deciding whether S is free is also NP-complete?
    # No. If "not free" is NP-complete, then "free" is co-NP-complete. A problem
    # cannot be both NP-complete and co-NP-complete unless NP = co-NP.
    answer_b = "No"

    # (c) If the problem of deciding whether S is free is in NP, and the problem of
    # deciding whether S is not free is NP-hard, does that imply the problem of
    # deciding whether S is free is NP-complete?
    # Yes. The given conditions imply NP = co-NP.
    # 1. "free" is in NP (given).
    # 2. "not free" is NP-hard (given).
    # From (1), "not free" is in co-NP.
    # From (2) and the conclusion that "not free" is in co-NP, it means an NP-hard problem
    # is in co-NP, which implies NP = co-NP.
    # Since "not free" is NP-hard, "free" is co-NP-hard.
    # If NP = co-NP, then co-NP-hard is the same as NP-hard.
    # So "free" is NP-hard.
    # Since "free" is in NP and is NP-hard, it is NP-complete.
    answer_c = "Yes"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_complexity_questions()