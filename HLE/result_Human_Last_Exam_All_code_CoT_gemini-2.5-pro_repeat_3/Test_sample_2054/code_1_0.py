def solve_complexity_questions():
    """
    Solves the three yes/no questions based on computational complexity theory.
    """
    # Part (a): If NOT-FREE is NP-hard, is FREE also NP-hard?
    # No. This would imply NP=co-NP, which is not proven.
    # If a problem P is NP-hard, its complement co-P is co-NP-hard.
    # Being co-NP-hard does not imply being NP-hard unless NP = co-NP.
    answer_a = "No"

    # Part (b): If NOT-FREE is NP-complete, is FREE also NP-complete?
    # No. If NOT-FREE is NP-complete, it is in NP and NP-hard.
    # Its complement, FREE, would then be in co-NP and co-NP-hard,
    # making it co-NP-complete. A co-NP-complete problem is not
    # NP-complete unless NP = co-NP.
    answer_b = "No"

    # Part (c): If FREE is in NP, and NOT-FREE is NP-hard, is FREE NP-complete?
    # Yes. The premises (FREE in NP implies NOT-FREE in co-NP; NOT-FREE is NP-hard)
    # together imply that NP = co-NP.
    # If NOT-FREE is NP-hard, then FREE is co-NP-hard.
    # If NP = co-NP, then co-NP-hard is the same as NP-hard.
    # So, FREE is NP-hard. Since it's also given to be in NP, it is NP-complete.
    answer_c = "Yes"

    final_answer = f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}"
    print(final_answer)

solve_complexity_questions()