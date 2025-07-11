def solve_enclitic_order():
    """
    Determines the correct order of a given list of Old Russian enclitics.
    """
    # This list represents the established canonical order of enclitics.
    # It is simplified for the purpose of this problem.
    # 1. Particles (бо, же)
    # 2. Conditional (бы)
    # 3. Pronouns (мя)
    # 4. Verbs (еси)
    canonical_order = [
        "бо", "же", "бы", "мя", "еси"
    ]

    # The list of enclitics to be sorted, as provided in the problem.
    given_enclitics = ["бы", "же", "еси", "бо", "мя"]

    # Sort the given list based on the index of each item in the canonical order.
    # The `key` function looks up the position of each enclitic in the master list.
    sorted_enclitics = sorted(given_enclitics, key=lambda clitic: canonical_order.index(clitic))

    print("If the enclitics were attached to the same word, their correct order would be:")
    # The problem asks to output each number in the final equation.
    # Here, 'number' is interpreted as each element in the sequence.
    print(" -> ".join(sorted_enclitics))

solve_enclitic_order()