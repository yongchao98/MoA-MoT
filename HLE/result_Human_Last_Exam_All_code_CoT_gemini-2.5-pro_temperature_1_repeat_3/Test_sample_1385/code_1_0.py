def solve_enclitics_order():
    """
    This function determines and prints the correct order of a given set of
    Old Russian enclitics based on established linguistic rules.
    """
    # The given enclitics: бы, же, еси, бо, мя
    # Based on linguistic rules, the order is:
    # 1. Particles (in the order: же, бы, бо)
    # 2. Pronominal clitics (мя)
    # 3. Verbal clitics (еси)
    ordered_enclitics = ["же", "бы", "бо", "мя", "еси"]
    
    print("The correct order for the enclitics is:")
    print(", ".join(ordered_enclitics))

solve_enclitics_order()
<<<же, бы, бо, мя, еси>>>