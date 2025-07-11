def solve_set_theory_problem():
    """
    Analyzes the relationship between an inaccessible cardinal and a measurable cardinal
    to determine the status of statement P within system S.
    """

    # The system S is ZFC + "There exists an inaccessible cardinal κ".
    # The statement P is "There exists a nontrivial elementary embedding j : V → M
    # with critical point κ".

    # Step 1: Translate the statement P into a more common term.
    # The existence of such an embedding 'j' is the definition of κ being a "measurable cardinal".
    # So, P is equivalent to the statement "κ is a measurable cardinal".
    explanation = """
1.  The Problem in Plain Terms:
    -   System S = ZFC + Axiom 'There exists an inaccessible cardinal κ'.
    -   Statement P = 'The cardinal κ is a measurable cardinal'.
    -   The question is: Can the system S prove P, disprove P, or is P independent of S?

2.  The Hierarchy of Large Cardinals:
    -   In set theory, large cardinal axioms are ranked by consistency strength.
    -   The axiom 'There exists a measurable cardinal' is known to be strictly stronger than 'There exists an inaccessible cardinal'.
    -   This means that if we assume ZFC + 'measurable' is consistent, it implies ZFC + 'inaccessible' is also consistent, but the reverse is not true.

3.  Conclusion:
    -   Because P (measurable) is a strictly stronger assertion than the axioms in S (inaccessible), S does not have the axiomatic power to prove P.
    -   Because it is consistent to have a measurable cardinal (which is also inaccessible), S cannot disprove P.
    -   When a statement can be neither proven nor disproven from a set of axioms, it is said to be 'independent' of them.
    -   Therefore, P is independent of S. This corresponds to choice C.
"""

    print(explanation)

    # The final answer is C.
    final_answer = "C"
    print(f"The correct choice is: {final_answer}")


solve_set_theory_problem()
<<<C>>>