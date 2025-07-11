def solve_enclitic_order():
    """
    Determines and prints the correct order of a given list of Old Russian enclitics.
    """
    # The list of enclitics provided by the user.
    enclitics_to_order = ['бы', 'же', 'еси', 'бо', 'мя']

    # Based on historical linguistics, the order is determined by the type of particle.
    # We can create a ranking map to define this order.
    # Lower numbers come first.
    # 1. же (conjunction)
    # 2. бы (conditional)
    # 3. бо (other particle)
    # 4. мя (pronoun)
    # 5. еси (verb)
    order_rules = {
        'же': 1,
        'бы': 2,
        'бо': 3,
        'мя': 4,
        'еси': 5
    }

    # Sort the list using the 'key' argument, which looks up the rank of each enclitic
    # in our 'order_rules' dictionary.
    sorted_enclitics = sorted(enclitics_to_order, key=lambda enclitic: order_rules[enclitic])

    # Print the final ordered list.
    # The prompt asks to output each "number", which in this context means each item in the sequence.
    print("The correct order for the enclitics is:")
    final_string = " → ".join(sorted_enclitics)
    print(final_string)

solve_enclitic_order()
