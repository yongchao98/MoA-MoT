def find_glider_species():
    """
    Identifies dragonfly species from a predefined list that are known
    for long-distance migration and a gliding lifestyle, an ecology
    associated with reduced pterostigmata.
    
    The identified species are:
    3) Macrodiplax balteata (Tribe Pantaliini)
    4) Pantala flavescens (Tribe Pantaliini - the 'Wandering Glider')
    10) Tholymis tillarga (Tribe Pantaliini)
    """
    
    # Indices of the species with the specified ecology.
    # These are 1-based indices as given in the problem description.
    glider_indices = [3, 4, 10]
    
    # To fulfill the requirement of outputting each number, we can show them
    # before joining.
    print("The indices of the species expected to have reduced pterostigmata are:")
    for index in glider_indices:
        print(f"Species index: {index}")
        
    # The final output should be the indices separated by a comma.
    # We convert the numbers to strings to use the join method.
    final_output_string = ",".join(map(str, glider_indices))
    
    print("\nFinal comma-separated list:")
    print(final_output_string)

find_glider_species()