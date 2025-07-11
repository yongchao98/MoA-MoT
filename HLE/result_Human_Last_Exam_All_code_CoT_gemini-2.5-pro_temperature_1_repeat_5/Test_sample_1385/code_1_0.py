def solve_enclitic_order():
    """
    Determines and prints the correct grammatical order of a given set of Old Russian enclitics.
    """

    # The list of enclitics provided by the user.
    # User input: бы, же, еси, бо, мя

    # The correct order is determined by Old Russian grammatical rules (Wackernagel's Law).
    # The hierarchy of enclitic classes is:
    # 1. Conjunctions/particles (же, бо)
    # 2. Pronominal clitics (мя)
    # 3. Conditional particle (бы)
    # 4. Verbal clitics (еси)
    # Within the first group, 'же' precedes 'бо'.

    correct_order = ["же", "бо", "мя", "бы", "еси"]

    print("The correct grammatical order for these enclitics is:")
    
    # We use join to create a formatted string from the list elements.
    # This fulfills the requirement to output each element of the final sequence.
    final_sequence = ", ".join(correct_order)
    print(final_sequence)

solve_enclitic_order()