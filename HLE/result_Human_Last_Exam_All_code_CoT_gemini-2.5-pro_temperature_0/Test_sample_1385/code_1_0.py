def solve_enclitic_order():
    """
    Determines and prints the correct grammatical order of a given list of
    Old Russian enclitics.
    """
    # The list of enclitics provided by the user.
    user_enclitics = ["бы", "же", "еси", "бо", "мя"]

    # The correct grammatical order of these enclitics, based on linguistic rules.
    # 1. же (conjunctive particle)
    # 2. бы (conditional particle)
    # 3. бо (causal particle)
    # 4. мя (enclitic pronoun)
    # 5. еси (enclitic verb form)
    correct_order = ["же", "бы", "бо", "мя", "еси"]

    # To sort the user's list, we use the index of each item in the
    # 'correct_order' list as the sorting key.
    # A dictionary mapping each enclitic to its correct position makes this clear.
    order_map = {enclitic: index for index, enclitic in enumerate(correct_order)}

    # Sort the user's list based on the position defined in the map.
    sorted_enclitics = sorted(user_enclitics, key=lambda e: order_map[e])

    # Print the final ordered list.
    print("The correct order for the enclitics is:")
    print(" -> ".join(sorted_enclitics))

solve_enclitic_order()
<<<же, бы, бо, мя, еси>>>