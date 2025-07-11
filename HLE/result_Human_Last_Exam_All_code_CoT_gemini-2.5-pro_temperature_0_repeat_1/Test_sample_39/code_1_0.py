def analyze_large_cardinals():
    """
    This function analyzes the relationship between inaccessible and measurable cardinals
    to determine the status of statement P within the formal system S.
    """

    # Step 1: Define the system and the statement
    system_S = "ZFC + 'There exists an inaccessible cardinal κ'"
    statement_P_description = "There exists a nontrivial elementary embedding j: V -> M with critical point κ"

    # Step 2: Interpret the statement in the context of large cardinal theory
    statement_P_equivalent = "κ is a measurable cardinal"

    print("Problem Analysis:")
    print(f"System S = {system_S}")
    print(f"Statement P asserts: '{statement_P_description}'")
    print(f"This is the definition of a measurable cardinal. So, P is equivalent to: '{statement_P_equivalent}'.")
    print("-" * 20)

    # Step 3: State the key theorem from set theory
    print("Key Concept: The Large Cardinal Hierarchy")
    print("The consistency strength of large cardinal axioms is well-ordered.")
    print("Fact: The consistency of 'there exists a measurable cardinal' is strictly greater than the consistency of 'there exists an inaccessible cardinal'.")
    print("-" * 20)

    # Step 4: Evaluate the options based on this fact
    print("Evaluating the Options:")

    # Option A: P is provable in S.
    print("A) Is P provable in S?")
    print("   No. If you could prove a stronger axiom (P) from a weaker set of axioms (S), they would have the same consistency strength. This contradicts the known hierarchy. So, A is false.")

    # Option B/D: The negation of P is provable in S / P leads to a contradiction.
    print("\nB/D) Is P disprovable in S or contradictory?")
    print("   No. This would imply that the existence of an inaccessible cardinal contradicts the existence of a measurable one. However, any model of ZFC with a measurable cardinal is also a model with an inaccessible cardinal (since every measurable is inaccessible). Thus, the theory 'S + P' is consistent if 'ZFC + P' is consistent. So, B and D are false.")

    # Option E: P is meaningless in S.
    print("\nE) Is P meaningless?")
    print("   No. The statement P is a well-formed and highly significant formula in the language of set theory. So, E is false.")

    # Option C: P is independent of S.
    print("\nC) Is P independent of S?")
    print("   Yes. Since P is neither provable nor disprovable from the axioms in S, it is by definition independent of S.")
    print("-" * 20)

    # Step 5: Final Conclusion
    final_answer = "C"
    print(f"Conclusion: The correct statement is that P is independent of S.")

analyze_large_cardinals()
<<<C>>>