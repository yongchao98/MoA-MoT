def solve_enclitic_order():
    """
    Determines and prints the correct order of a given set of Old Russian enclitics.
    """
    # The user's list of enclitics.
    enclitics_to_order = ["бы", "же", "еси", "бо", "мя"]

    # The canonical order of these enclitics, based on historical linguistic rules.
    # 1. же (emphatic/conjunctive particle)
    # 2. бы (conditional particle)
    # 3. мя (accusative personal pronoun)
    # 4. еси (enclitic verb form, 'you are')
    # 5. бо (causal particle, 'for')
    correct_order = ["же", "бы", "мя", "еси", "бо"]

    # Sort the user's list based on the position of each item in the correct_order list.
    # The `key` argument of the `sorted` function looks up the index of each clitic
    # in the `correct_order` list to determine its sorting position.
    sorted_enclitics = sorted(enclitics_to_order, key=lambda clitic: correct_order.index(clitic))

    # Print the final, correctly ordered list.
    print("The correct order for the given enclitics is:")
    print(", ".join(sorted_enclitics))

solve_enclitic_order()