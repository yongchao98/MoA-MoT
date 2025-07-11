def solve_complexity_questions():
    """
    This function provides answers to the complexity theory questions
    based on the relationships between complexity classes. The reasoning
    is detailed in the comments.
    """

    # Let's define the problems for clarity:
    # P1: The problem "deciding whether S is not free".
    # P2: The problem "deciding whether S is free".
    # P2 is the complement of P1 (P2 = co-P1).

    # (a) If the problem of deciding whether S is not free is NP-hard,
    # does that imply the problem of deciding whether S is free is also NP-hard?

    # Analysis for (a):
    # Given: P1 is NP-hard.
    # Question: Is P2 NP-hard?
    # By definition, if a problem is NP-hard, its complement is co-NP-hard.
    # So, if P1 is NP-hard, then P2 is co-NP-hard.
    # A problem being co-NP-hard does not imply it is NP-hard unless NP = co-NP,
    # which is a major open problem in computer science. Without this assumption,
    # the implication does not hold.
    answer_a = "No"

    # (b) If the problem of deciding whether S is not free is NP-complete,
    # does that imply the problem of deciding whether S is free is also NP-complete?

    # Analysis for (b):
    # Given: P1 is NP-complete. This means P1 is in NP and P1 is NP-hard.
    # Question: Is P2 NP-complete?
    # For P2 to be NP-complete, it must be in NP and be NP-hard.
    # 1. Is P2 in NP? Since P1 is in NP, its complement P2 is in co-NP. P2 would
    #    be in NP only if NP = co-NP.
    # 2. Is P2 NP-hard? As per (a), since P1 is NP-hard, P2 is co-NP-hard. It would
    #    be NP-hard only if NP = co-NP.
    # Since the implication relies on the unproven assumption NP = co-NP, it does not hold.
    # If P1 is NP-complete, P2 is co-NP-complete.
    answer_b = "No"

    # (c) If the problem of deciding whether S is free is in NP, and the problem
    # of deciding whether S is not free is NP-hard, does that imply the problem
    # of deciding whether S is free is NP-complete?

    # Analysis for (c):
    # Given:
    # 1. P2 is in NP.
    # 2. P1 is NP-hard.
    # Question: Is P2 NP-complete?
    # For P2 to be NP-complete, it must be (i) in NP and (ii) NP-hard.
    # (i) Condition 1 gives us that P2 is in NP.
    # (ii) Let's check if P2 is NP-hard. From condition 2 (P1 is NP-hard),
    #      we know its complement, P2, is co-NP-hard.
    # We now have that P2 is in NP and P2 is co-NP-hard.
    # A known theorem in complexity theory states that if a co-NP-hard problem is
    # also in NP, it implies that NP = co-NP.
    # If NP = co-NP, then the classes NP-hard and co-NP-hard are identical.
    # Since P2 is co-NP-hard, it must also be NP-hard.
    # Both conditions for P2 being NP-complete are satisfied. Therefore, the implication is true.
    answer_c = "Yes"

    # Print the final answer in the required format.
    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_complexity_questions()