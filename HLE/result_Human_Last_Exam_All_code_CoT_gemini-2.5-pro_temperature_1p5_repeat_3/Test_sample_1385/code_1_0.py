def solve_enclitic_order():
    """
    This function determines and prints the correct order of a given set of
    Old Russian enclitics based on established linguistic rules.
    """

    # The established order of Old Russian enclitics based on historical linguistics.
    # The hierarchy is generally:
    # 1. Conjunctive/Emphatic particles (e.g., же, бо)
    # 2. Auxiliary verb clitics (e.g., еси)
    # 3. Pronominal clitics (e.g., мя)
    # 4. Conditional particle (e.g., бы)

    # The enclitics arranged in their correct final order
    ordered_enclitics = ["же", "бо", "еси", "мя", "бы"]

    print("The correct order for the enclitics is:")
    
    # Print each element of the final ordered sequence
    final_sequence = " -> ".join(ordered_enclitics)
    print(final_sequence)

solve_enclitic_order()