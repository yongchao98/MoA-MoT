def find_glider_dragonflies():
    """
    Identifies dragonflies expected to have reduced pterostigmata based on their ecology.

    This function is based on the biological principle that long-distance migratory
    dragonflies that utilize a gliding flight style tend to have reduced pterostigmata
    to minimize aerodynamic drag.
    """
    
    # The full list of taxa provided by the user.
    all_taxa = [
        "Didymops transversa",      # 1
        "Urothemis edwarsi",        # 2
        "Macrodiplax balteata",     # 3
        "Pantala flavescens",       # 4
        "Orthetrum cancellatum",    # 5
        "Libelulla quadrimaculata", # 6
        "Libelulla pulchela",       # 7
        "Sympetrum corruptum",      # 8
        "Celithemis elisa",         # 9
        "Tholymis tillarga"         # 10
    ]

    # A set of species known for their gliding/migratory ecology.
    # These are the species expected to have reduced pterostigmata.
    glider_species = {
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Tholymis tillarga"
    }

    # A list to store the 1-based indices of the matching species.
    result_indices = []

    # Iterate through the list of all taxa with their index.
    for index, taxon in enumerate(all_taxa):
        if taxon in glider_species:
            # Add the 1-based index to our results.
            result_indices.append(str(index + 1))
            
    # Check if we found any matching species and format the output.
    if result_indices:
        final_answer = ",".join(result_indices)
    else:
        final_answer = "none"
        
    print(final_answer)

find_glider_dragonflies()