def order_old_russian_enclitics():
    """
    This script determines the correct grammatical order for a set of
    Old Russian enclitics and prints the result.
    """
    # The list of enclitics provided in the problem.
    enclitics_to_sort = ["бы", "же", "еси", "бо", "мя"]

    # The correct order is based on established rules of Old Russian linguistics.
    # We create a dictionary to map each enclitic to its rank in the sequence.
    # Lower numbers come first.
    # 1. же (particle)
    # 2. бо (particle, follows же)
    # 3. бы (conditional particle)
    # 4. еси (auxiliary verb)
    # 5. мя (pronoun)
    correct_order_map = {
        "же": 1,
        "бо": 2,
        "бы": 3,
        "еси": 4,
        "мя": 5
    }

    # We sort the original list using the ranks from our map as the sorting key.
    sorted_enclitics = sorted(enclitics_to_sort, key=lambda word: correct_order_map[word])

    # Print the final ordered list. Each word is an element in the final sequence.
    print("The correct order would be:")
    # We join with ' - ' to simulate how they would attach to a preceding word.
    print("WORD-" + " - ".join(sorted_enclitics))
    
    # As requested, here is each "item" from the final sequence:
    print("\nThe individual enclitics in order are:")
    for item in sorted_enclitics:
        print(item)


order_old_russian_enclitics()
<<<же, бо, бы, еси, мя>>>