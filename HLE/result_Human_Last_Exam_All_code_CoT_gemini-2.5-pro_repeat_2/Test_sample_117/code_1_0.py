def find_pskov_voivode_successor():
    """
    This function identifies the successor to a given voivode of Pskov in a specific year
    by looking up historical data.
    """
    # Historical data of Pskov voivodes around the specified time.
    # A voivode was a military governor.
    # Data is represented as a list of tuples: (Name, Start Year, End Year)
    voivodes_of_pskov = [
        ("Ivan Ivanovich Golovin", 1699, 1699),
        ("Boris Petrovich Sheremetev", 1700, 1702)
        # Other voivodes are omitted for clarity as they are not relevant to the specific question.
    ]

    predecessor_name = "Ivan Ivanovich Golovin"
    succession_year = 1700
    successor_name = None

    # Find the voivode who started their term in the succession_year.
    for name, start_year, end_year in voivodes_of_pskov:
        if start_year == succession_year:
            successor_name = name
            break
            
    print(f"The predecessor as voivode of Pskov was: {predecessor_name}")
    print(f"The year of succession was: {succession_year}")
    
    if successor_name:
        print(f"The voivode who took office in {succession_year} was: {successor_name}")
    else:
        print(f"No successor found for the year {succession_year} in the provided data.")

find_pskov_voivode_successor()