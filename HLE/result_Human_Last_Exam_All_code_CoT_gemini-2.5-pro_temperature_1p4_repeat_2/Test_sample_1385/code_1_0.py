def print_enclitic_order():
    """
    This function defines the correct order of the given Old Russian enclitics
    and prints them in sequence.
    """
    # The established grammatical order of the enclitics
    # 1. бо (particle)
    # 2. же (particle)
    # 3. бы (conditional particle)
    # 4. мя (pronominal clitic)
    # 5. еси (verbal clitic)
    enclitics_in_order = ['бо', 'же', 'бы', 'мя', 'еси']

    # Print the final ordered sequence
    print("The correct order of the enclitics is:")
    # We output each "number" (item) in the final "equation" (sequence) as requested.
    print(" -> ".join(enclitics_in_order))

print_enclitic_order()