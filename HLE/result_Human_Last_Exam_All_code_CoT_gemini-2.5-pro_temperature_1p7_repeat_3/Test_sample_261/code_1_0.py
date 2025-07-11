def find_archimandrite():
    """
    Finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period based on historical data.
    """
    # Historical data of archimandrites around the specified time.
    # Data is in the format: (Name, Start Year, End Year)
    archimandrite_list = [
        ("Amvrosiy", 1729, 1730),
        ("Markell", 1730, 1731),
        ("Innokentiy", 1731, 1736)
    ]

    target_start_year = 1730
    target_end_year = 1731
    
    found_name = None

    for name, start, end in archimandrite_list:
        if start == target_start_year and end == target_end_year:
            found_name = name
            break
            
    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from the year {target_start_year} to {target_end_year} was {found_name}.")
    else:
        print(f"Could not find the archimandrite for the period {target_start_year}-{target_end_year} in the provided data.")

find_archimandrite()