def solve_set_theory_problem():
    """
    Analyzes the relationship between inaccessible and measurable cardinals in ZFC.

    The problem sets up a system S and a statement P:
    - System S: ZFC + "There exists an inaccessible cardinal κ".
    - Statement P: "There exists a nontrivial elementary embedding j : V → M such that the critical point of j is κ".

    This analysis follows these steps:
    1.  Statement P is the definition of "κ is a measurable cardinal".
    2.  The existence of a measurable cardinal is a much stronger large cardinal axiom than the existence of an inaccessible cardinal. A theory cannot prove a statement that is strictly stronger than itself in consistency strength. Therefore, S cannot prove P.
    3.  Conversely, S cannot disprove P. The consistency of S with the negation of P (¬P) can be shown by considering Gödel's constructible universe L. If S is consistent, then a model for S + ¬P exists.
    4.  Since S can neither prove P nor disprove P, the statement P is independent of the system S.
    """
    explanation = solve_set_theory_problem.__doc__
    print("--- Reasoning ---")
    print(explanation)
    print("--- Conclusion ---")
    print("The statement P is independent of the system S.")
    print("The correct choice is C.")

solve_set_theory_problem()

# The final answer is determined by the logical independence of the two large cardinal axioms.
# The existence of a measurable cardinal (Statement P) is a stronger axiom than the
# existence of an inaccessible cardinal (assumed in System S). Therefore, S can neither
# prove nor disprove P.
# The final answer is C.
print("<<<C>>>")