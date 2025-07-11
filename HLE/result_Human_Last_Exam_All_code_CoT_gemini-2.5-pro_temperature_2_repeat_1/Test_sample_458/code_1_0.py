def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The primary ecological trait associated with reduced pterostigmata is a
    migratory, gliding flight style, as this life history selects for
    maximum flight efficiency over the anti-flutter properties provided by a
    larger pterostigma in flapping species.

    The species selected are known long-distance migrants and gliders:
    - 3: Macrodiplax balteata (Marl Pennant)
    - 4: Pantala flavescens (Globe Skimmer)
    - 8: Sympetrum corruptum (Variegated Meadowhawk)
    - 10: Tholymis tillarga (Coral-tailed Cloudwing)
    """

    # Indices of the species with migratory/gliding ecologies
    candidate_indices = [3, 4, 8, 10]

    # Convert the list of integers to a comma-separated string for the output
    result = ",".join(map(str, candidate_indices))

    print(result)

find_species_with_reduced_pterostigmata()