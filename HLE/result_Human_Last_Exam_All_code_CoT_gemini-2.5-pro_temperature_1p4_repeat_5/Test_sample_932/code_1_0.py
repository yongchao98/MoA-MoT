def solve_beat_sheet_task():
    """
    Determines which insect immatures are unlikely to be collected via a beat-sheet method.

    The beat-sheet method works by dislodging insects living on the exterior of plants.
    Therefore, we need to identify tribes whose immatures live in other habitats.
    """
    tribes = [
        {"index": 1, "name": "Apis", "immature_habitat": "Internal nest (hive with wax cells)"},
        {"index": 2, "name": "Melipotini", "immature_habitat": "External on plants (caterpillars on foliage)"},
        {"index": 3, "name": "Eupholini", "immature_habitat": "Internal borer (weevil larvae in plant tissue)"},
        {"index": 4, "name": "Acritini", "immature_habitat": "In decaying matter (predatory larvae in dung/carrion)"},
        {"index": 5, "name": "Oxyptilini", "immature_habitat": "External on plants (caterpillars on foliage/flowers)"},
        {"index": 6, "name": "Dictyophorini", "immature_habitat": "External on plants (nymphs on stems/leaves)"},
        {"index": 7, "name": "Acanthocerini", "immature_habitat": "In decaying wood (larvae in rotting logs)"}
    ]

    unlikely_indices = []

    for tribe in tribes:
        # The beat-sheet method collects insects from the exterior of plants.
        # If the habitat is NOT 'External on plants', they are unlikely to be collected.
        if "External on plants" not in tribe["immature_habitat"]:
            unlikely_indices.append(tribe["index"])

    # Sort the indices in ascending order as requested.
    unlikely_indices.sort()
    
    # Format the output as a comma-separated string.
    # The str(i) converts each integer index to a string.
    # ", ".join(...) combines them into a single string.
    final_answer = ", ".join(map(str, unlikely_indices))
    
    print(final_answer)

solve_beat_sheet_task()