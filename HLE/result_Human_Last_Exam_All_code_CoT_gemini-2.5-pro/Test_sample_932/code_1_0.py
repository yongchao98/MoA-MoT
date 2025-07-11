def find_unlikely_tribes():
    """
    Identifies insect tribes whose immatures are unlikely to be collected
    using a beat-sheet method based on their habitat.
    """
    # A dictionary representing the tribes and the typical habitat of their immatures.
    # Habitats suitable for beat-sheet collection are 'foliage' or 'stems'.
    tribes_info = {
        1: {'name': 'Apis', 'habitat': 'in hive'},
        2: {'name': 'Melipotini', 'habitat': 'foliage'},
        3: {'name': 'Eupholini', 'habitat': 'foliage'},
        4: {'name': 'Acritini', 'habitat': 'decaying matter'},
        5: {'name': 'Oxyptilini', 'habitat': 'foliage'},
        6: {'name': 'Dictyophorini', 'habitat': 'stems'},
        7: {'name': 'Acanthocerini', 'habitat': 'decaying wood'}
    }

    unlikely_indices = []
    
    print("Analyzing the habitats of the immature stages of each tribe:")
    for index, data in tribes_info.items():
        # The beat-sheet method is effective for insects on external plant parts.
        is_on_foliage = data['habitat'] in ['foliage', 'stems']
        
        if is_on_foliage:
            print(f"{index}) {data['name']}: Immatures are found on foliage/stems. Likely to be collected.")
        else:
            print(f"{index}) {data['name']}: Immatures are found {data['habitat']}. Unlikely to be collected.")
            unlikely_indices.append(index)

    # Sort the indices in ascending order
    unlikely_indices.sort()

    print("\nThe indices of the tribes unlikely to be collected are:")
    # Print each number separated by a comma, as requested.
    # The * operator unpacks the list elements as separate arguments to the print function.
    print(*unlikely_indices, sep=", ")

find_unlikely_tribes()