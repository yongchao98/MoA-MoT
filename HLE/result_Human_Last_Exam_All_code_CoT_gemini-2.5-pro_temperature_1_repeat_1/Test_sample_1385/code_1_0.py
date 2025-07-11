def solve_enclitic_order():
    """
    This function determines and prints the correct grammatical order of a given
    list of Old Russian enclitics.
    """

    # The list of enclitics to be ordered: бы, же, еси, бо, мя
    
    # Based on linguistic rules, the correct order is:
    # 1. Particles (же, бы, бо)
    # 2. Enclitic verb (еси)
    # 3. Enclitic pronoun (мя)
    correctly_ordered_enclitics = ['же', 'бы', 'бо', 'еси', 'мя']

    print("The correct order for the enclitics, if they were attached to the same word, is:")
    
    # The prompt requests to output each word in a final "equation".
    # We will join the words with ' + ' to create this format.
    final_equation = " + ".join(correctly_ordered_enclitics)
    
    print(final_equation)

solve_enclitic_order()