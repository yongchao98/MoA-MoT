def solve_enclitic_order():
    """
    This function determines and prints the correct order of a given
    list of Old Russian enclitics based on established linguistic rules.
    """

    # The list of enclitics to be ordered.
    # The user provided: бы, же, еси, бо, мя
    
    # The correct order is determined by linguistic hierarchy:
    # 1. Conjunctional particle: бо
    # 2. Focusing particle: же
    # 3. Conditional particle: бы
    # 4. Pronominal enclitic: мя
    # 5. Verbal enclitic: еси
    correct_order = ["бо", "же", "бы", "мя", "еси"]

    # Combine the enclitics into a single string to represent them
    # attached to a word. This fulfills the request to output each
    # element in the final sequence.
    final_sequence = "-".join(correct_order)

    print("The correct order of the enclitics is:")
    print(final_sequence)

solve_enclitic_order()