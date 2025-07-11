def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a list that are expected to have
    reduced pterostigmata based on their ecology.

    The primary ecological factor is a specialization for long-distance
    migratory flight using gliding and soaring, which creates an evolutionary
    pressure to reduce wing weight and drag.
    """

    # The full list of taxa provided
    all_taxa = [
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

    # Species known for their "glider" ecology (long-distance migration and soaring)
    # which is associated with reduced pterostigmata.
    target_species = [
        "Macrodiplax balteata",
        "Pantala flavescens",
        "Tholymis tillarga"
    ]

    # A list to hold the 1-based indices of the matching species
    result_indices = []

    # Iterate through the list of all taxa to find the indices of our target species
    for index, species_name in enumerate(all_taxa):
        if species_name in target_species:
            # Add the 1-based index (index + 1) to the results
            result_indices.append(str(index + 1))

    # Join the found indices with a comma and print the result
    final_output = ",".join(result_indices)
    print(final_output)

find_species_with_reduced_pterostigmata()