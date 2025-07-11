def solve_enclitic_order():
    """
    Determines and prints the fixed order of a given set of Old Russian enclitics.
    """
    # The order of enclitics in Old Russian followed a strict template.
    # The hierarchy, based on linguistic analysis, is generally as follows:
    # 1. Causal/Conjunctive particles (like 'бо')
    # 2. Emphatic/Contrastive particles (like 'же')
    # 3. Conditional particle ('бы')
    # 4. Pronominal clitics (like the accusative 'мя')
    # 5. Verbal clitics (auxiliary verbs like 'еси')

    # The given enclitics are: 'бы', 'же', 'еси', 'бо', 'мя'
    # Applying the hierarchy gives the following order:
    ordered_enclitics = ["бо", "же", "бы", "мя", "еси"]

    # We will print this order in the form of an equation, as requested.
    output_string = " + ".join(ordered_enclitics)
    
    print("Based on the grammatical rules of Old Russian, the correct order for these enclitics is:")
    print(output_string)

solve_enclitic_order()