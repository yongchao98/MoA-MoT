def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The primary ecological factor for reduced pterostigmata is a lifestyle
    of long-distance migration and extensive gliding, where the aerodynamic
    drag of the pterostigma outweighs its anti-flutter benefits.
    """

    # Species data provided by the user.
    # Key: index, Value: [Name, Ecology Notes]
    species_info = {
        1: ["Didymops transversa", "Stream percher, not a specialized glider."],
        2: ["Urothemis edwarsi", "Percher near still water, not a specialized glider."],
        3: ["Macrodiplax balteata", "Strong long-distance migrant and glider."],
        4: ["Pantala flavescens", "Archetypal long-distance migrant and high-altitude glider."],
        5: ["Orthetrum cancellatum", "Territorial percher, not a specialized glider."],
        6: ["Libelulla quadrimaculata", "Migratory, but not a specialized glider; has prominent pterostigmata."],
        7: ["Libelulla pulchela", "Territorial percher, not a specialized glider."],
        8: ["Sympetrum corruptum", "Migratory, but not primarily a glider."],
        9: ["Celithemis elisa", "Territorial percher, not a specialized glider."],
        10: ["Tholymis tillarga", "Long-distance migrant and crepuscular glider."]
    }

    species_indices_with_reduced_pterostigmata = []
    
    # We identify the species known as long-distance gliders.
    # Based on established entomological knowledge, these are Pantala,
    # Tholymis, and Macrodiplax.
    target_species_keywords = ["glider", "gliding"]

    for index, data in species_info.items():
        name = data[0]
        ecology_note = data[1]
        
        # Check if the ecology note contains keywords related to gliding migration.
        if any(keyword in ecology_note for keyword in target_species_keywords):
            species_indices_with_reduced_pterostigmata.append(index)

    # Sort the indices for consistent output
    species_indices_with_reduced_pterostigmata.sort()

    # Format the result as a comma-separated string
    result_string = ",".join(map(str, species_indices_with_reduced_pterostigmata))

    print(result_string)

find_species_with_reduced_pterostigmata()