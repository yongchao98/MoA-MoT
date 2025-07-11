def solve_enclitic_order():
    """
    Determines and prints the correct order of specified Old Russian enclitics.
    """
    # The established order of Old Russian enclitics based on linguistic rules.
    # 1. же (emphatic/conjunctive particle)
    # 2. бы (conditional particle)
    # 3. мя (pronominal clitic)
    # 4. еси (verbal clitic)
    # 5. бо (causal particle)
    correct_order = ["же", "бы", "мя", "еси", "бо"]

    # Format the output as a sequence.
    # The prompt asks to "output each number in the final equation", which we interpret
    # as outputting each element in the sequence.
    output_string = " + ".join(correct_order)

    print("The correct order of the enclitics is:")
    print(output_string)

solve_enclitic_order()