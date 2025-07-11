def solve_beat_sheet_task():
    """
    This function determines which insect tribes' immatures are unlikely to
    be collected using a beat-sheet method and prints their indices.
    """

    # Dictionary of tribes with reasoning for collection likelihood.
    # The beat-sheet method primarily collects insects living on the exterior of plants.
    tribes = {
        1: ("Apis", "Unlikely: Immatures are in hives, not on foliage."),
        2: ("Melipotini", "Likely: Immatures are foliage-feeding caterpillars."),
        3: ("Eupholini", "Unlikely: Immatures are internal wood-boring larvae."),
        4: ("Acritini", "Unlikely: Immatures are predators in dung, carrion, or nests."),
        5: ("Oxyptilini", "Unlikely: Immatures are often internal stem or flower borers."),
        6: ("Dictyophorini", "Likely: Immatures (nymphs) live on plant stems and leaves."),
        7: ("Acanthocerini", "Unlikely: Immatures are soil-dwelling or wood-boring grubs.")
    }

    # Identify the indices of tribes unlikely to be collected.
    unlikely_indices = []
    for index, (name, reason) in tribes.items():
        if "Unlikely" in reason:
            unlikely_indices.append(index)

    # Sort the indices in ascending order.
    unlikely_indices.sort()

    # Format the output string as comma-separated values.
    # The str(num) part ensures each number in the final equation is outputted.
    output_string = ", ".join([str(num) for num in unlikely_indices])
    
    # Print the final result.
    print(output_string)

solve_beat_sheet_task()