def solve_beat_sheet_puzzle():
    """
    Determines which insect tribes' immatures are unlikely to be collected
    using a beat-sheet method based on their life history.

    A beat-sheet is effective for insects living externally on plant foliage.
    """

    tribes_data = {
        # Index: (Tribe Name, Immature Habitat, Collectible via Beat-Sheet)
        1: ("Apis", "In a nest/hive", False),
        2: ("Melipotini", "External on foliage (caterpillar)", True),
        3: ("Eupholini", "Internal plant borer (grub)", False),
        4: ("Acritini", "In dung/carrion/litter (predator)", False),
        5: ("Oxyptilini", "External on foliage (caterpillar)", True),
        6: ("Dictyophorini", "External on foliage (nymph)", True),
        7: ("Acanthocerini", "In decaying wood/litter (grub)", False)
    }

    unlikely_to_be_collected_indices = []
    for index, data in tribes_data.items():
        is_collectible = data[2]
        if not is_collectible:
            unlikely_to_be_collected_indices.append(index)

    # Sort the indices in ascending order
    unlikely_to_be_collected_indices.sort()

    # Format the result as a string with commas
    # The map(str, ...) converts each integer index to a string
    result_string = ", ".join(map(str, unlikely_to_be_collected_indices))

    print(result_string)

solve_beat_sheet_puzzle()