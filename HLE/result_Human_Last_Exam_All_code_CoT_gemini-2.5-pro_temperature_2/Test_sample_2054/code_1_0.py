def solve_complexity_questions():
    """
    Analyzes and solves the three complexity theory questions.
    """

    # Let P be the decision problem "S is free".
    # Let co-P be its complement, "S is not free".

    # (a) If the problem of deciding whether S is not free is NP-hard,
    # does that imply the problem of deciding whether S is free is also NP-hard?
    #
    # Reasoning:
    # We are given that co-P is NP-hard. The question asks if P must be NP-hard.
    # By definition, if a problem L is NP-hard, its complement co-L is co-NP-hard.
    # So, if co-P is NP-hard, then P is co-NP-hard.
    # The question is whether being co-NP-hard implies being NP-hard.
    # This is only true if NP = co-NP, which is a major unresolved question in computer science
    # and is widely believed to be false. Unless we assume NP = co-NP, we cannot conclude this implication holds.
    # For example, assuming NP != co-NP, the TAUTOLOGY problem is co-NP-hard but not NP-hard.
    # Therefore, the implication is not generally valid.
    answer_a = "No"

    # (b) If the problem of deciding whether S is not free is NP-complete,
    # does that imply the problem of deciding whether S is free is also NP-complete?
    #
    # Reasoning:
    # We are given that co-P is NP-complete. This means co-P is in NP and is NP-hard.
    # The question asks if P is NP-complete.
    # If co-P is NP-complete, its complement P is, by definition, co-NP-complete.
    # A problem is co-NP-complete if it is in co-NP and is co-NP-hard.
    # For P to be NP-complete, it must be in NP and be NP-hard.
    # If a problem were both NP-complete and co-NP-complete, it would imply that NP = co-NP.
    # As this is not known to be true, we cannot conclude that a co-NP-complete problem is also NP-complete.
    # For instance, if SAT is NP-complete, its complement TAUTOLOGY is co-NP-complete, but not known to be NP-complete.
    # Therefore, the implication is not valid.
    answer_b = "No"

    # (c) If the problem of deciding whether S is free is in NP, and the problem
    # of deciding whether S is not free is NP-hard, does that imply the problem
    # of deciding whether S is free is NP-complete?
    #
    # Reasoning:
    # We are given two premises:
    # 1. P is in NP.
    # 2. co-P is NP-hard.
    # We need to determine if this implies P is NP-complete. To do this, we need to show
    # that P is NP-hard (since we already know P is in NP from premise 1).
    #
    # From premise 1 (P is in NP), its complement, co-P, must be in co-NP.
    # So, we have a problem, co-P, that is both in co-NP (from premise 1) and is NP-hard (from premise 2).
    # A fundamental theorem of complexity theory states that if an NP-hard problem is also in co-NP, then NP = co-NP.
    # Thus, the premises together imply that NP = co-NP.
    #
    # Now, let's show P is NP-hard.
    # From premise 2, co-P is NP-hard. This means any NP problem reduces to it. For example, SAT <=p co-P.
    # By taking the complements, this reduction implies co-SAT <=p P.
    # Since co-SAT is co-NP-complete, this reduction proves that P is co-NP-hard.
    # Since the premises imply NP = co-NP, the classes of NP-hard and co-NP-hard problems are identical.
    # Therefore, P being co-NP-hard means it is also NP-hard.
    #
    # Since P is in NP (premise 1) and P is NP-hard (as we have shown), P is NP-complete by definition.
    # Thus, the implication is valid.
    answer_c = "Yes"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}")

# Execute the function to print the answers.
solve_complexity_questions()