def solve_complexity_questions():
    """
    Solves the three given theoretical questions about computational complexity classes.

    The questions are about the relationship between a problem and its complement
    in terms of NP-hardness and NP-completeness.

    Let P_free be the problem "is the semigroup S free?".
    Let P_not_free be its complement "is the semigroup S not free?".

    (a) If P_not_free is NP-hard, is P_free also NP-hard?
    The class of NP-hard problems is not known to be closed under complement.
    If P_not_free is NP-hard, then P_free is co-NP-hard. These are different
    unless NP = co-NP, which is not known. Thus, the implication does not hold.
    Answer: No.

    (b) If P_not_free is NP-complete, is P_free also NP-complete?
    If P_not_free is NP-complete, it lies in NP. Its complement, P_free, must then
    lie in co-NP. For P_free to be NP-complete, it must be in NP. This would only
    happen if NP = co-NP. If NP != co-NP, the complement of an NP-complete problem
    cannot be NP-complete. Thus, the implication does not hold.
    Answer: No.

    (c) If P_free is in NP, and P_not_free is NP-hard, does that imply
    P_free is NP-complete?
    The premises are that P_free is in NP and its complement is NP-hard.
    These premises are strong enough to imply that NP = co-NP. In a world where
    NP = co-NP, the property of being NP-hard is equivalent to being co-NP-hard.
    The fact that P_not_free is NP-hard implies P_free is co-NP-hard, and thus
    NP-hard. Since P_free is in NP (given) and is also NP-hard, it is NP-complete
    by definition. The implication holds.
    Answer: Yes.
    """
    answer_a = "No"
    answer_b = "No"
    answer_c = "Yes"

    print(f"(a) {answer_a}; (b) {answer_b}; (c) {answer_c}.")

solve_complexity_questions()