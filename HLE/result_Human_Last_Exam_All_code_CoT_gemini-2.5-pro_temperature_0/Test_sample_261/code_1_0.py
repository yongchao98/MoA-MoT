def find_archimandrite():
    """
    This function finds the archimandrite of the Pskov-Caves Monastery
    for a specific period by searching through a predefined list of data.
    """
    # Data representing some of the archimandrites and their service years.
    # Format: (Name, Start Year, End Year)
    monastery_leaders = [
        ('Korniliy', 1529, 1570),
        ('Markell', 1727, 1731),
        ('Veniamin', 1731, 1732),
        ('Amvrosiy', 1740, 1744),
        ('Innokentiy', 1850, 1868)
    ]

    # The period we are interested in.
    start_year_query = 1730
    end_year_query = 1731

    print(f"Searching for the archimandrite who served from {start_year_query} to {end_year_query}...")

    found_leader = None
    for name, start_year, end_year in monastery_leaders:
        # Check if the leader's tenure covers the entire query period.
        # A leader is a match if their service started on or before the query start
        # and ended on or after the query end.
        if start_year <= start_year_query and end_year >= end_year_query:
            found_leader = name
            print(f"\nFound a match: {name}")
            print(f"Service period: {start_year} - {end_year}")
            print(f"This period includes the queried years: {start_year_query} - {end_year_query}")
            break

    if not found_leader:
        print("\nNo archimandrite found in the dataset for the specified period.")

find_archimandrite()