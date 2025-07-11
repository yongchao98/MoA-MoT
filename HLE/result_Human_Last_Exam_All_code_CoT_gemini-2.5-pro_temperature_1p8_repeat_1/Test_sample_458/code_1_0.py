def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.
    """
    species_list = [
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

    # Species known for long-distance migration and extensive gliding flight,
    # an ecology that favors the reduction of the pterostigma to minimize drag.
    migratory_gliders = [
        "Pantala flavescens",
        "Tholymis tillarga"
    ]

    print("Identifying species with reduced pterostigmata based on ecology (long-distance gliding/migration):")
    
    # Find the 1-based indices of the matching species
    indices = []
    for i, species in enumerate(species_list):
        if species in migratory_gliders:
            # The list is 0-indexed, but the question uses 1-based indices.
            indices.append(str(i + 1))
            print(f"- Found: {species} (Index {i+1})")

    # Format the final result string
    if not indices:
        result = "none"
    else:
        result = ",".join(indices)

    print("\nFinal list of indices for species expected to have reduced pterostigmata:")
    print(result)

find_species_with_reduced_pterostigmata()