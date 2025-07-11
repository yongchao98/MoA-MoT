def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected to have 
    reduced pterostigmata based on their ecology.

    The primary ecological factor is long-distance migratory behavior, characteristic
    of species often called "gliders". These species, such as Pantala flavescens,
    Macrodiplax balteata, and Tholymis tillarga, undertake extensive migrations,
    often across oceans. This lifestyle is correlated with a reduction in the size
    of the pterostigma, a feature that typically aids in flight stability at high speeds 
    but may be less critical or energetically costly for prolonged gliding flight.

    The function returns the 1-based indices of these species.
    """
    
    species_list = [
        "Didymops transversa",
        "Urothemis edwarsi",
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Orthetrum cancellatum",
        "Libelulla quadrimaculata",
        "Libelulla pulchela",
        "Sympetrum corruptum",
        "Celithemis elisa",
        "Tholymis tillarga"
    ]

    # Indices (1-based) of species known for long-distance migration
    # and associated reduced pterostigmata.
    # 3: Macrodiplax balteata (Marl Pennant)
    # 4: Pantala flavescens (Wandering Glider)
    # 10: Tholymis tillarga (Coral-tailed Cloudwing)
    migratory_glider_indices = [3, 4, 10]

    # Format the indices as a comma-separated string for the final output
    # Using list comprehension and str.join for an efficient and clean implementation.
    output_string = ",".join([str(index) for index in migratory_glider_indices])

    print(output_string)

find_species_with_reduced_pterostigmata()
<<<3,4,10>>>