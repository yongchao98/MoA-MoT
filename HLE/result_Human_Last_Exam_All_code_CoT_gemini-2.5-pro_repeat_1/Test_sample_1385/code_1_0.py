def solve_enclitic_order():
    """
    Determines and prints the correct grammatical order for a given
    set of Old Russian enclitics.
    """
    # The list of enclitics to be ordered, as provided by the user.
    enclitics_to_order = ['бы', 'же', 'еси', 'бо', 'мя']

    # The established correct order of these enclitics based on historical linguistics.
    # This list acts as a key for sorting. The order is:
    # 1. же (emphatic particle)
    # 2. бы (conditional particle)
    # 3. бо (causal particle)
    # 4. мя (pronominal clitic - accusative)
    # 5. еси (verbal clitic - 2nd p. sg. of 'to be')
    correct_order_key = ['же', 'бы', 'бо', 'мя', 'еси']

    # We sort the initial list. The 'key' argument uses the index of each
    # word in our 'correct_order_key' list to determine its sort position.
    # For example, 'же' is at index 0 in the key, so it will come first.
    # 'бы' is at index 1, so it will come second, and so on.
    sorted_enclitics = sorted(enclitics_to_order, key=correct_order_key.index)

    # Print the final, correctly ordered list of enclitics.
    print("If attached to the same word, the enclitics would appear in the following order:")
    
    # We join the elements of the list into a single string for clean output.
    final_order_string = ", ".join(sorted_enclitics)
    print(final_order_string)

solve_enclitic_order()

<<<же, бы, бо, мя, еси>>>