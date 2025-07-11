def solve_enclitic_order():
    """
    Determines and prints the correct grammatical order for a given list of
    Old Russian enclitics.
    """
    # The list of enclitics provided by the user.
    enclitics_to_sort = ['бы', 'же', 'еси', 'бо', 'мя']

    # The correct order of enclitics is based on established linguistic rules.
    # We can represent this order by assigning a rank to each enclitic.
    # Lower rank means it appears earlier in the sequence.
    # 1. Particles: бо (for), же (indeed), бы (would)
    # 2. Pronouns: мя (me)
    # 3. Verbs: еси (are)
    # The internal order of the particles is also fixed: бо -> же -> бы.
    correct_order_map = {
        'бо': 1,
        'же': 2,
        'бы': 3,
        'мя': 4,
        'еси': 5
    }

    # We sort the given list using the ranks defined in our map.
    # The 'key' argument uses a lambda function to look up the rank of each word.
    sorted_enclitics = sorted(enclitics_to_sort, key=lambda word: correct_order_map[word])

    # The prompt asks to output each "number" in the final "equation".
    # We will interpret this as printing each word in the final sequence.
    print("The initial list of enclitics is:", enclitics_to_sort)
    print("The correct order of these enclitics is:")
    
    # We print the final ordered list as a sequence.
    final_sequence_string = ", ".join(sorted_enclitics)
    print(final_sequence_string)

solve_enclitic_order()
