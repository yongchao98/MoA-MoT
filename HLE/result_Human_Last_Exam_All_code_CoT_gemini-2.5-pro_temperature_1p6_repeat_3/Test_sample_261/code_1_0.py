def find_archimandrite():
    """
    Finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period based on a predefined list.
    """
    # A list of archimandrites and their years of service around the target period.
    # Data is based on historical records of the monastery's leadership.
    archimandrites_data = [
        {"name": "Gennadiy II", "start": 1727, "end": 1729},
        {"name": "Markell", "start": 1730, "end": 1731},
        {"name": "Gedeon", "start": 1731, "end": 1736}
    ]

    target_start_year = 1730
    target_end_year = 1731
    found_archimandrite = None

    for person in archimandrites_data:
        if person["start"] == target_start_year and person["end"] == target_end_year:
            found_archimandrite = person
            break

    if found_archimandrite:
        name = found_archimandrite["name"]
        start_year = found_archimandrite["start"]
        end_year = found_archimandrite["end"]
        print(f"The archimandrite of the Pskov-Caves Monastery from {start_year} to {end_year} was {name}.")
    else:
        print(f"Could not find the archimandrite for the period {target_start_year}-{target_end_year} in the provided data.")

find_archimandrite()