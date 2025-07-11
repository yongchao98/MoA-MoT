def solve_beat_sheet_task():
    """
    Analyzes insect tribes to determine which are unlikely to be collected
    as immatures using a beat-sheet method.

    The beat-sheet method is effective for collecting insects that live freely on
    the external surfaces of plants (leaves, stems). This script identifies
    tribes whose immatures live in other habitats.
    """

    # Data on the lifestyle of the immature stage of each tribe.
    # 'collectable' is True if immatures live on plant foliage.
    tribes_data = {
        1: {'name': 'Apis', 'reason': 'Immatures live in a protected hive/nest, not on plants.', 'collectable': False},
        2: {'name': 'Melipotini', 'reason': 'Immatures (caterpillars) are foliage feeders.', 'collectable': True},
        3: {'name': 'Eupholini', 'reason': 'Immatures are internal borers within plant tissues.', 'collectable': False},
        4: {'name': 'Acritini', 'reason': 'Immatures are predators in decaying matter/leaf litter, not on live foliage.', 'collectable': False},
        5: {'name': 'Oxyptilini', 'reason': 'Immatures (caterpillars) feed on leaves and flowers.', 'collectable': True},
        6: {'name': 'Dictyophorini', 'reason': 'Immatures (nymphs) live and feed externally on plants.', 'collectable': True},
        7: {'name': 'Acanthocerini', 'reason': 'Immatures live within dung balls, often buried underground.', 'collectable': False}
    }

    unlikely_indices = []
    
    # Iterate through the tribes and identify those unlikely to be collected.
    for index, data in tribes_data.items():
        if not data['collectable']:
            unlikely_indices.append(index)
    
    # The indices are already in order, but sorting ensures correctness.
    unlikely_indices.sort()
    
    # Format the final result as a comma-separated string.
    final_answer = ", ".join(map(str, unlikely_indices))

    print("The beat-sheet method collects insects living freely on plant foliage.")
    print("Based on the lifestyles of their immatures, the following tribes are unlikely to be collected:")
    for index in unlikely_indices:
        print(f" - Tribe {index} ({tribes_data[index]['name']}): {tribes_data[index]['reason']}")

    print("\nFinal list of indices:")
    print(final_answer)

solve_beat_sheet_task()