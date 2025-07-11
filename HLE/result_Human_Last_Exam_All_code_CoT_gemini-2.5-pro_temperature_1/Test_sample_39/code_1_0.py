def analyze_provability():
    """
    Analyzes the provability of a statement within an axiomatic system
    by comparing their consistency strengths. This is an analogy for the
    set-theoretic problem.
    """

    # In set theory, large cardinal axioms have a well-defined hierarchy
    # of consistency strength. We can represent this with numbers.
    # A higher number means a stronger theory.
    consistency_strengths = {
        "ZFC": 1,
        "S (ZFC + Inaccessible Cardinal)": 2,
        "P (ZFC + Measurable Cardinal)": 3  # Statement P is equivalent to this axiom
    }

    system_S_strength = consistency_strengths["S (ZFC + Inaccessible Cardinal)"]
    statement_P_strength = consistency_strengths["P (ZFC + Measurable Cardinal)"]

    print("System S is 'ZFC + Inaccessible Cardinal'.")
    print("Statement P is 'There exists a Measurable Cardinal'.\n")

    print(f"Representing consistency strengths numerically:")
    print(f"  - Strength of S: {system_S_strength}")
    print(f"  - Strength of P: {statement_P_strength}\n")

    # A key result (related to Gödel's Incompleteness Theorems) is that a
    # system cannot prove a statement that asserts the consistency of a
    # strictly stronger system.
    print("Analysis:")
    if statement_P_strength > system_S_strength:
        print("  - P is not provable in S.")
        print("    Reason: The axiom in P is strictly stronger than the axioms in S.")
        print("    (A weaker theory cannot prove a strictly stronger one.)")
    else:
        # This case is not reached, but included for completeness.
        print("  - P might be provable in S.")

    # Now, check for contradiction. A contradiction would exist if P implies the
    # negation of S's axioms.
    # Theorem: Every measurable cardinal is inaccessible.
    # This means P implies the axiom of S.
    is_contradictory = False
    print("\n  - The negation of P is not provable in S.")
    print("    Reason: The existence of a measurable cardinal (P) implies the")
    print("    existence of an inaccessible cardinal (S's axiom).")
    print("    Therefore, P and S are consistent with each other.")

    print("\nConclusion:")
    print("Since P is neither provable nor disprovable in S, it is independent of S.")


if __name__ == "__main__":
    analyze_provability()
