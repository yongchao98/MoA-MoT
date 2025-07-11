def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species expected to have reduced pterostigmata based on their ecology.

    The primary ecological driver for a reduced pterostigma is a lifestyle
    of long-distance, high-altitude gliding migration, where reducing drag
    and weight is paramount.
    """
    species_data = {
        1: {"name": "Didymops transversa", "ecology": ["stream-dweller", "percher"]},
        2: {"name": "Urothemis edwarsi", "ecology": ["pond-dweller", "local migrant"]},
        3: {"name": "Macrodiplax balteata", "ecology": ["long-distance migrant", "coastal", "glider"]},
        4: {"name": "Pantala flavescens", "ecology": ["global migrant", "ocean-crosser", "high-altitude glider"]},
        5: {"name": "Orthetrum cancellatum", "ecology": ["pond-dweller", "percher", "local disperser"]},
        6: {"name": "Libelulla quadrimaculata", "ecology": ["strong flier", "large-scale migrant"]},
        7: {"name": "Libelulla pulchela", "ecology": ["territorial", "percher"]},
        8: {"name": "Sympetrum corruptum", "ecology": ["migrant", "pond-dweller"]},
        9: {"name": "Celithemis elisa", "ecology": ["pond-dweller", "percher"]},
        10: {"name": "Tholymis tillarga", "ecology": ["crepuscular", "migrant", "glider"]},
    }

    # Define the key ecological traits associated with reduced pterostigmata
    target_traits = ["glider", "high-altitude glider", "ocean-crosser"]

    # A list to hold the indices of the identified species
    result_indices = []

    for index, data in species_data.items():
        # Check if any of the species' ecological traits match our target traits
        if any(trait in data["ecology"] for trait in target_traits):
            result_indices.append(str(index))

    # Join the indices with a comma and print the result
    output = ",".join(result_indices)
    print(output)

find_species_with_reduced_pterostigmata()