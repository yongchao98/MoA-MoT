def find_unlikely_tribes():
    """
    Identifies insect tribes whose immatures are unlikely to be collected
    using a beat-sheet method based on their life history.
    """

    # Data on the typical habitat of immatures for each tribe.
    # 'external_foliage' means they are likely to be collected by a beat-sheet.
    # Other habitats mean they are unlikely to be collected.
    tribes = {
        1: {'name': 'Apis', 'immature_habitat': 'hive'},
        2: {'name': 'Melipotini', 'immature_habitat': 'external_foliage'},
        3: {'name': 'Eupholini', 'immature_habitat': 'internal_borer_or_soil'},
        4: {'name': 'Acritini', 'immature_habitat': 'dung_or_carrion'},
        5: {'name': 'Oxyptilini', 'immature_habitat': 'internal_borer'},
        6: {'name': 'Dictyophorini', 'immature_habitat': 'external_foliage'},
        7: {'name': 'Acanthocerini', 'immature_habitat': 'soil_or_rotting_wood'}
    }

    unlikely_indices = []
    for index, data in tribes.items():
        # A beat-sheet collects insects from external foliage.
        # If the habitat is not 'external_foliage', it's unlikely to be collected.
        if data['immature_habitat'] != 'external_foliage':
            unlikely_indices.append(index)

    # The list is already in ascending order because we iterated through the sorted dictionary keys.
    # Convert each integer index to a string for joining.
    unlikely_indices_str = [str(i) for i in unlikely_indices]
    
    # Print the final list of numbers as a comma-separated string.
    print(", ".join(unlikely_indices_str))

find_unlikely_tribes()