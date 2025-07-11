def solve_beat_sheet_puzzle():
    """
    Determines which insect immatures are unlikely to be collected via a beat-sheet method.

    The function encodes the life history of each tribe and filters for those
    whose immatures live in concealed habitats rather than freely on foliage.
    """
    # 1. Define the tribes and the typical habitat of their immatures.
    #    'external' means on leaves/stems, making them collectible by beat-sheet.
    #    'internal_or_concealed' means they are inside nests, soil, wood, or dung,
    #    making them uncollectible by this method.
    tribes_info = {
        1: {'name': 'Apis', 'habitat': 'internal_or_concealed'},        # In a hive
        2: {'name': 'Melipotini', 'habitat': 'external'},               # Caterpillars on leaves
        3: {'name': 'Eupholini', 'habitat': 'internal_or_concealed'},   # Larvae in soil or stems
        4: {'name': 'Acritini', 'habitat': 'internal_or_concealed'},    # Larvae in dung/carrion
        5: {'name': 'Oxyptilini', 'habitat': 'external'},               # Caterpillars on leaves/flowers
        6: {'name': 'Dictyophorini', 'habitat': 'external'},            # Nymphs on stems/leaves
        7: {'name': 'Acanthocerini', 'habitat': 'internal_or_concealed'}# Grubs in soil/wood
    }

    # 2. Identify the indices of tribes with concealed immatures.
    unlikely_to_be_collected_indices = []
    for index, info in tribes_info.items():
        if info['habitat'] == 'internal_or_concealed':
            unlikely_to_be_collected_indices.append(index)

    # 3. Sort the indices in ascending order (they are already sorted but it's good practice).
    unlikely_to_be_collected_indices.sort()

    # 4. Print the result as a comma-separated string.
    result_string = ','.join(map(str, unlikely_to_be_collected_indices))
    print(result_string)

solve_beat_sheet_puzzle()