def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a list that are known to have
    reduced pterostigmata based on their ecology.

    The key ecological trait is being a high-altitude "glider" that relies
    on soaring for long-distance migration. This flight style favors
    adaptations to reduce wing mass and drag, which includes the reduction
    of the pterostigma.
    """
    
    # The list of taxa provided by the user
    taxa = [
        "Didymops transversa",       # 1
        "Urothemis edwarsi",         # 2
        "Macrodiplax balteata",      # 3
        "Pantala flavescens",        # 4
        "Orthetrum cancellatum",     # 5
        "Libelulla quadrimaculata",  # 6
        "Libelulla pulchela",        # 7
        "Sympetrum corruptum",       # 8
        "Celithemis elisa",          # 9
        "Tholymis tillarga"          # 10
    ]

    # Species known for a "glider" ecology and reduced pterostigmata
    glider_species_with_reduced_pterostigmata = {
        "Pantala flavescens",  # The Globe Skimmer, a classic example
        "Tholymis tillarga"   # The Coral-tailed Cloudwing, also a known glider
    }

    result_indices = []
    # Iterate through the list and find the 1-based indices of matching species
    for index, species_name in enumerate(taxa):
        if species_name in glider_species_with_reduced_pterostigmata:
            result_indices.append(str(index + 1))

    if result_indices:
        print(",".join(result_indices))
    else:
        print("none")

# Execute the function to get the answer
find_species_with_reduced_pterostigmata()