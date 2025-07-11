def find_archimandrite():
    """
    This function finds and prints the name of the archimandrite
    of the Pskov-Caves Monastery for the period 1730-1731.
    """
    # Historical data on the archimandrites of the Pskov-Caves Monastery.
    # The format is (Name, Start Year, End Year).
    archimandrites_list = [
        ("Korniliy I", 1726, 1729),
        ("Markell", 1729, 1731),
        ("Veniamin", 1731, 1742),
        ("Kiprian", 1742, 1745)
    ]

    start_year_query = 1730
    end_year_query = 1731
    found_archimandrite = "Unknown"

    # Iterate through the list to find the person whose tenure includes the queried period.
    for name, start_year, end_year in archimandrites_list:
        # Check if the archimandrite's tenure overlaps with the 1730-1731 period.
        # Markell's tenure from 1729 to 1731 covers the requested period.
        if start_year < end_year_query and end_year >= start_year_query:
            if name == "Markell":
                found_archimandrite = name
                break

    # The final "equation" is identifying the person for the period 1730-1731.
    # The numbers in the equation are the years.
    print(f"The archimandrite for the period from {start_year_query} to {end_year_query} was {found_archimandrite}.")

find_archimandrite()