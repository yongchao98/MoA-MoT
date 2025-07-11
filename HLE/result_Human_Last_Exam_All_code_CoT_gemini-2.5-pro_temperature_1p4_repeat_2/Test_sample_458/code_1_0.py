def find_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that likely have reduced pterostigmata
    based on their flight ecology.
    """
    
    species_list = [
        "1) Didymops transversa (Stream Cruiser)",
        "2) Urothemis edwarsi (Black-tailed Skimmer)",
        "3) Macrodiplax balteata (Marl Pennant)",
        "4) Pantala flavescens (Wandering Glider)",
        "5) Orthetrum cancellatum (Black-tailed Skimmer)",
        "6) Libelulla quadrimaculata (Four-spotted Chaser)",
        "7) Libelulla pulchela (Twelve-spotted Skimmer)",
        "8) Sympetrum corruptum (Variegated Meadowhawk)",
        "9) Celithemis elisa (Calico Pennant)",
        "10) Tholymis tillarga (Coral-tailed Cloudwing)"
    ]

    # Indices of species expected to have reduced pterostigmata.
    # We start with an empty list and add the indices based on ecological reasoning.
    reduced_pterostigmata_indices = []

    # Analysis:
    # 4) Pantala flavescens: A world-traveling "glider". This species is a supreme
    #    soarer, part of the Pantalini tribe. These dragonflies have very broad
    #    hindwings for gliding and relatively small pterostigmata for their size.
    reduced_pterostigmata_indices.append(4)

    # 9) Celithemis elisa: A "pennant". Pennants are known for their relatively
    #    weak, fluttering flight, perching frequently. This flight style does not
    #    necessitate the powerful wing stabilization a large pterostigma provides.
    reduced_pterostigmata_indices.append(9)

    # 10) Tholymis tillarga: Another member of the Pantalini tribe (the "gliders"),
    #     like Pantala. It is also a migratory glider with broad hindwings and
    #     a comparatively small pterostigma adapted for that lifestyle.
    reduced_pterostigmata_indices.append(10)

    # Print the explanation for the chosen species
    print("Species expected to have reduced pterostigmata based on their ecology:")
    print("4: Pantala flavescens is a specialized, long-distance glider.")
    print("9: Celithemis elisa is a weak, fluttering flier (a 'pennant').")
    print("10: Tholymis tillarga is also a specialized glider, closely related to Pantala.")
    
    # Format the final answer as a comma-separated string
    answer = ",".join(map(str, sorted(reduced_pterostigmata_indices)))
    
    print("\nThe indices of the species are:")
    print(answer)

    # Final answer in the specified format
    print(f"<<<{answer}>>>")

find_reduced_pterostigmata()