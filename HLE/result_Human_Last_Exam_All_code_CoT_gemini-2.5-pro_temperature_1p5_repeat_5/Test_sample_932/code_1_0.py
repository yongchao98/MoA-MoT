def solve_beat_sheet_task():
    """
    Determines which insect tribes' immatures are unlikely to be collected
    using a beat-sheet method based on their typical habitat.
    """

    # Data: Tribe name and the typical habitat of their immatures.
    # 'External Foliage' means they are likely to be collected by a beat-sheet.
    # 'Other' includes internal borers, hive-dwellers, soil-dwellers, etc.
    tribes_data = {
        1: {'name': 'Apis', 'habitat': 'Other'},          # Live in hives
        2: {'name': 'Melipotini', 'habitat': 'External Foliage'}, # External leaf-feeders
        3: {'name': 'Eupholini', 'habitat': 'Other'},         # Internal borers
        4: {'name': 'Acritini', 'habitat': 'Other'},          # Predators in dung/carrion
        5: {'name': 'Oxyptilini', 'habitat': 'Other'},        # Internal borers in buds/stems
        6: {'name': 'Dictyophorini', 'habitat': 'External Foliage'}, # External sap-feeders
        7: {'name': 'Acanthocerini', 'habitat': 'Other'}      # Live in soil/rotting wood
    }

    unlikely_to_be_collected_indices = []
    for index, data in tribes_data.items():
        if data['habitat'] == 'Other':
            unlikely_to_be_collected_indices.append(index)

    # Sort the indices in ascending order
    unlikely_to_be_collected_indices.sort()

    # Print the result as a comma-separated string
    # The map(str, ...) converts each integer in the list to a string
    # The ",".join(...) concatenates the string elements with a comma
    print(",".join(map(str, unlikely_to_be_collected_indices)))

solve_beat_sheet_task()