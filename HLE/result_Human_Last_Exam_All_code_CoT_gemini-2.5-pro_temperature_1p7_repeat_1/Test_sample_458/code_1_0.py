def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected to have
    reduced pterostigmata based on their ecology.

    This determination is based on the principle that reduced pterostigmata are an
    adaptation for efficient long-distance gliding flight. Species known for
    migration and a "glider" lifestyle are the primary candidates.
    """

    species_list = {
        1: "Didymops transversa",      # Cruiser, not a primary long-distance glider
        2: "Urothemis edwarsi",        # Percher
        3: "Macrodiplax balteata",     # Migratory "wanderer", strong candidate
        4: "Pantala flavescens",       # Quintessential "glider" and trans-oceanic migrant
        5: "Orthetrum cancellatum",    # Percher
        6: "Libelulla quadrimaculata", # Migratory, but more of a chaser/percher style
        7: "Libelulla pulchela",       # Percher
        8: "Sympetrum corruptum",      # Migratory, but generally a meadowhawk/percher
        9: "Celithemis elisa",         # Percher
        10: "Tholymis tillarga"        # Migratory "cloudwing", high-flyer/glider
    }

    # Rationale for selection:
    # We are looking for species whose ecology is dominated by long-distance,
    # gliding flight. The pterostigma adds weight that helps prevent wing flutter
    # but can be less optimal for maximum gliding efficiency.
    # Therefore, highly migratory "gliders" and "wanderers" are expected to
    # exhibit reduced pterostigmata.
    #
    # - Macrodiplax balteata (#3) is a known migratory "wanderer".
    # - Pantala flavescens (#4) is the world's most famous "wandering glider".
    # - Tholymis tillarga (#10) is a migratory "cloudwing" known for high-altitude flight.

    expected_species_indices = [3, 4, 10]
    
    # Sort indices for consistent output
    expected_species_indices.sort()

    # Convert indices to strings and join with a comma
    result = ",".join(map(str, expected_species_indices))
    
    print("Based on their ecology as long-distance migrants and gliders, the following species are expected to have reduced pterostigmata:")
    print("The indices of these species are:")
    print(result)

find_species_with_reduced_pterostigmata()