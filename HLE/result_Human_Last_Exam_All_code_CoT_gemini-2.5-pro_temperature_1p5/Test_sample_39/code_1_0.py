import textwrap

def explain_set_theory_problem():
    """
    Analyzes the relationship between an inaccessible cardinal and the existence
    of an elementary embedding to determine the status of statement P.
    """

    # Step 1: Define the system and the statement
    print("Step 1: Analyzing the problem components.")
    print("System S = ZFC + 'There exists an inaccessible cardinal κ'.")
    print("Statement P = 'There exists a nontrivial elementary embedding j: V → M with critical point κ'.")
    print("-" * 30)

    # Step 2: Translate statement P into a known large cardinal property
    print("Step 2: Translating statement P.")
    explanation_p = """
    The statement P is equivalent to a standard large cardinal axiom. A cardinal κ is the critical point of a nontrivial elementary embedding j: V → M if and only if κ is what is known as a 'measurable cardinal'. So, P is equivalent to the statement 'The cardinal κ is a measurable cardinal'.
    """
    print(textwrap.dedent(explanation_p).strip())
    print("-" * 30)

    # Step 3: Compare the consistency strength of the axioms
    print("Step 3: Comparing the strength of inaccessible and measurable cardinals.")
    explanation_strength = """
    In set theory, large cardinal axioms are ordered by their 'consistency strength'. The existence of a measurable cardinal is a much stronger assumption than the existence of an inaccessible cardinal.

    - Every measurable cardinal is inaccessible.
    - However, it is not provable from the existence of an inaccessible cardinal that it (or any other cardinal) must be measurable.
    """
    print(textwrap.dedent(explanation_strength).strip())
    print("-" * 30)

    # Step 4: Determine the status of P within S
    print("Step 4: Determining the status of P in system S.")
    explanation_conclusion = """
    - Since 'measurable' is a stronger property than 'inaccessible', we cannot prove P within system S. Proving P would be like proving a stronger theory from a weaker one.
    - We also cannot disprove P within system S. The axioms of S are consistent with P (assuming ZFC + 'there exists a measurable cardinal' is itself consistent).
    - Therefore, P is neither provable nor disprovable from the axioms of S.
    """
    print(textwrap.dedent(explanation_conclusion).strip())
    print("-" * 30)

    # Step 5: Final Answer
    print("Conclusion: P is independent of S.")
    print("\nThis corresponds to answer choice C.")


if __name__ == "__main__":
    explain_set_theory_problem()