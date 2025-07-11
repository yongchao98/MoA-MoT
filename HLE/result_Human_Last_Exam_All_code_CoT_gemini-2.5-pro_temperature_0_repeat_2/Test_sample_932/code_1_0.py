def solve_beat_sheet_task():
    """
    Determines which insect tribes' immatures are unlikely to be collected
    by a beat-sheet method based on their life habits.

    The beat-sheet method collects insects living externally on foliage.
    This function identifies tribes whose immatures are internal borers,
    live in hives/nests, or in decaying matter.
    """

    # Data: Tribe index, name, and a boolean indicating if their immatures
    # live on foliage (True) or not (False).
    tribes_data = {
        1: {"name": "Apis", "on_foliage": False},           # Larvae are in hives
        2: {"name": "Melipotini", "on_foliage": True},      # Caterpillars on leaves
        3: {"name": "Eupholini", "on_foliage": False},      # Larvae are internal borers
        4: {"name": "Acritini", "on_foliage": False},       # Larvae in decaying matter
        5: {"name": "Oxyptilini", "on_foliage": True},      # Caterpillars on leaves/flowers
        6: {"name": "Dictyophorini", "on_foliage": True},   # Nymphs on stems/leaves
        7: {"name": "Acanthocerini", "on_foliage": False}   # Larvae in rotting wood/nests
    }

    unlikely_indices = []
    for index, data in tribes_data.items():
        if not data["on_foliage"]:
            unlikely_indices.append(index)

    # Sort the indices in ascending order as requested
    unlikely_indices.sort()

    # Format the output string
    result_string = ", ".join(map(str, unlikely_indices))

    print(f"The indices of the tribes unlikely to be collected are: {result_string}")

solve_beat_sheet_task()
<<<1, 3, 4, 7>>>