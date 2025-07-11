def find_species_with_reduced_pterostigmata():
    """
    Analyzes a list of dragonfly species to identify which is expected to have 
    a reduced pterostigma based on its flight ecology.

    The pterostigma prevents wing flutter at high speeds. Therefore:
    - Strong fliers and gliders need a large pterostigma.
    - Weak, fluttering fliers have less need for one, and it may be reduced.
    """
    
    species_list = [
        "Didymops transversa",      # 1: Strong flier (cruiser)
        "Urothemis edwarsi",        # 2: Strong flier (basker)
        "Macrodiplax balteata",     # 3: Strong migrant
        "Pantala flavescens",       # 4: Premier long-distance glider
        "Orthetrum cancellatum",    # 5: Strong flier (skimmer)
        "Libelulla quadrimaculata", # 6: Strong migrant (chaser)
        "Libelulla pulchela",       # 7: Strong flier (skimmer)
        "Sympetrum corruptum",      # 8: Strong migrant
        "Celithemis elisa",         # 9: Weak, fluttering flier (pennant)
        "Tholymis tillarga"         # 10: Strong migrant (glider)
    ]

    # The genus "Celithemis" (pennants) is noted for its weak, fluttering flight style.
    # This makes it the most likely candidate for reduced pterostigmata.
    target_genus = "Celithemis"
    
    reduced_pterostigmata_indices = []
    
    for i, species_name in enumerate(species_list):
        if species_name.startswith(target_genus):
            # Indices are 1-based in the user's list
            reduced_pterostigmata_indices.append(str(i + 1))

    if reduced_pterostigmata_indices:
        result = ",".join(reduced_pterostigmata_indices)
    else:
        result = "none"
        
    print(result)

find_species_with_reduced_pterostigmata()