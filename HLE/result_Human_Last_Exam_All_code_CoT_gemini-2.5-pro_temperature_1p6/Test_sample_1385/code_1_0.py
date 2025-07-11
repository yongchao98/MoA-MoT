def solve_enclitic_order():
    """
    Determines and prints the correct order of a given set of Old Russian enclitics.
    """
    # The given list of enclitics in an arbitrary order.
    given_enclitics = ["бы", "же", "еси", "бо", "мя"]

    # The canonical order of enclitics based on linguistic rules.
    # 1. бы (conditional)
    # 2. мя (pronominal)
    # 3. же (emphatic)
    # 4. бо (emphatic, follows же)
    # 5. еси (verbal)
    correct_order = ["бы", "мя", "же", "бо", "еси"]
    
    # We can "sort" the given list by iterating through the correct order
    # and picking the elements that are in our given list.
    final_ordered_list = []
    for particle in correct_order:
        if particle in given_enclitics:
            final_ordered_list.append(particle)

    # Print the final result as a string.
    print(f"The correct order for the enclitics is:")
    print(", ".join(final_ordered_list))

solve_enclitic_order()