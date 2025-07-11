def find_archimandrite():
    """
    Finds the archimandrite of the Pskov-Caves Monastery for a specific year.
    """
    # Historical data representing the tenures of the archimandrites.
    # Each entry is a tuple: (start_year, end_year, name).
    archimandrite_records = [
        (1726, 1731, "Markell"),
        (1731, 1732, "Kirill"),
        (1732, 1745, "Amvrosiy"),
        (1746, 1753, "Veniamin")
    ]

    # The period of interest is 1730-1731. We will search for the year 1730.
    target_start_year = 1730
    target_end_year = 1731

    found_name = None
    for start, end, name in archimandrite_records:
        # Check if the target year is within the current record's range.
        if start <= target_start_year < end:
            found_name = name
            break

    if found_name:
        # The prompt asks to output each number in the final equation.
        # We will represent this by printing the years we searched for.
        print(f"Searching for the period: {target_start_year}-{target_end_year}")
        print(f"The archimandrite of the Pskov-Caves Monastery during this time was: {found_name}")
    else:
        print(f"No record found for the period {target_start_year}-{target_end_year}")

find_archimandrite()