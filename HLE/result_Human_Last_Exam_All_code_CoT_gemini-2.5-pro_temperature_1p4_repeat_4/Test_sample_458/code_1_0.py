def find_species_with_reduced_pterostigmata():
    """
    Identifies dragonfly species from a predefined list that are expected
    to have reduced pterostigmata based on their ecology.

    The key ecological trait associated with reduced pterostigmata is being a
    highly adapted, long-distance migrant and glider. This code filters for
    species known for this behavior.
    """
    
    # A list of dictionaries representing the taxa and their relevant ecological traits.
    # The 'ecology' key notes if the species is a well-known long-distance migrant
    # whose flight style is associated with reduced pterostigmata.
    taxa_data = [
        {'index': 1, 'name': 'Didymops transversa', 'ecology': 'stream cruiser, non-migrant'},
        {'index': 2, 'name': 'Urothemis edwarsi', 'ecology': 'local migrant, not a glider'},
        {'index': 3, 'name': 'Macrodiplax balteata', 'ecology': 'long-distance migrant'},
        {'index': 4, 'name': 'Pantala flavescens', 'ecology': 'long-distance migrant'},
        {'index': 5, 'name': 'Orthetrum cancellatum', 'ecology': 'local colonizer'},
        {'index': 6, 'name': 'Libelulla quadrimaculata', 'ecology': 'irruptive migrant, prominent pterostigmata'},
        {'index': 7, 'name': 'Libelulla pulchela', 'ecology': 'percher, non-migrant'},
        {'index': 8, 'name': 'Sympetrum corruptum', 'ecology': 'long-distance migrant'},
        {'index': 9, 'name': 'Celithemis elisa', 'ecology': 'percher, non-migrant'},
        {'index': 10, 'name': 'Tholymis tillarga', 'ecology': 'long-distance migrant'}
    ]

    migrant_indices = []
    for species in taxa_data:
        # We are looking for species known as strong, long-distance migrants/gliders.
        if species['ecology'] == 'long-distance migrant':
            migrant_indices.append(str(species['index']))
    
    # Print the final result as a comma-separated string
    print(",".join(migrant_indices))

find_species_with_reduced_pterostigmata()