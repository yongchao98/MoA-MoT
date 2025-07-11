def solve_enclitic_order():
    """
    This function determines and prints the correct order of a given set of
    Old Russian enclitics based on established linguistic rules.
    """
    # The given enclitics: бы, же, еси, бо, мя
    
    # Based on linguistic hierarchy, the correct order is:
    # 1. бо (conjunction)
    # 2. же (emphatic particle)
    # 3. бы (conditional particle)
    # 4. еси (verbal clitic)
    # 5. мя (pronominal clitic)
    
    ordered_enclitics = ["бо", "же", "бы", "еси", "мя"]
    
    # Print the final ordered sequence.
    print("The correct order of the enclitics is:")
    print(" -> ".join(ordered_enclitics))

solve_enclitic_order()
<<<бо -> же -> бы -> еси -> мя>>>