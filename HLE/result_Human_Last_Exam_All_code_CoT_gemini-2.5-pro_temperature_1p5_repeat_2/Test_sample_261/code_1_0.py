def find_archimandrite():
    """
    Finds the archimandrite of the Pskov-Caves Monastery for the years 1730-1731.
    """
    # Historical data of the archimandrites and their years of service around the target period.
    # Format: (Name, Start Year, End Year)
    archimandrites_data = [
        ("Varlaam", 1727, 1729),
        ("Markell", 1730, 1731),
        ("Veniamin", 1731, 1738)
    ]

    target_start_year = 1730
    target_end_year = 1731
    found_name = None

    for name, start_year, end_year in archimandrites_data:
        # Check if the archimandrite's term matches the target period.
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break
            
    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_name}")
    else:
        print(f"No archimandrite found for the period {target_start_year}-{target_end_year} in the provided data.")

find_archimandrite()