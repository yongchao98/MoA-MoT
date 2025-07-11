def find_gliding_dragonflies():
    """
    Identifies dragonfly species from a given list that are expected to have
    reduced pterostigmata based on their gliding ecology.
    """

    # The pterostigma adds weight to the wing tip, preventing flutter during
    # fast, continuous flapping. Species that are specialized gliders or soarers
    # do not rely on this type of flight and thus often have reduced pterostigmata.
    species_data = {
        1: ("Didymops transversa", "Cruiser: A strong, patrolling flier that flies continuously."),
        2: ("Urothemis edwarsi", "Percher: Not a specialized long-distance glider."),
        3: ("Macrodiplax balteata", "Migrant: A strong flier over long distances, but not primarily a glider."),
        4: ("Pantala flavescens", "Wandering Glider: A famous trans-oceanic migrant known for its extensive use of gliding and soaring on air currents."),
        5: ("Orthetrum cancellatum", "Skimmer: A percher and territorial flier, not a specialized glider."),
        6: ("Libelulla quadrimaculata", "Chaser: A powerful and active flier."),
        7: ("Libelulla pulchela", "Skimmer: A territorial percher with strong, rapid flights."),
        8: ("Sympetrum corruptum", "Migrant: A powerful migratory flier."),
        9: ("Celithemis elisa", "Pennant: A small percher with short, fluttery flights."),
        10: ("Tholymis tillarga", "Cloudwing: A migratory species that flies high and is known to soar and glide, especially during migration.")
    }

    print("Analyzing dragonfly ecology to find species with reduced pterostigmata...")
    print("Principle: Gliding and soaring species require less wing-flutter damping and are expected to have reduced pterostigmata.\n")

    reduced_pterostigmata_indices = []
    for index, (name, ecology) in species_data.items():
        # Search for keywords indicating a gliding or soaring flight style.
        if "glider" in ecology.lower() or "soaring" in ecology.lower() or "soars" in ecology.lower():
            reduced_pterostigmata_indices.append(index)
            print(f"- Found candidate: Index {index}, {name}. Ecology: {ecology}")

    # Sort indices for consistent output
    reduced_pterostigmata_indices.sort()

    # Format the result as a comma-separated string
    result_string = ",".join(map(str, reduced_pterostigmata_indices))

    if not result_string:
        result_string = "none"

    print(f"\nFinal Answer: The indices of species expected to have reduced pterostigmata are {result_string}.")

find_gliding_dragonflies()
<<<4,10>>>