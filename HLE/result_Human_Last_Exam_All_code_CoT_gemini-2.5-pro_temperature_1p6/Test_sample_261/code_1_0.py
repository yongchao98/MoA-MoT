def find_archimandrite():
    """
    This function retrieves the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period from a predefined data set.
    """
    # A small dataset of the monastery's leaders based on historical records.
    # The keys are the periods of service and the values are the names.
    archimandrites_data = {
        "1725-1729": "Varlaam",
        "1730-1731": "Markell",
        "1732-1743": "Amvrosiy",
        "1743-1746": "Gennadiy"
    }

    # The period we are querying.
    start_year = 1730
    end_year = 1731
    query_period = f"{start_year}-{end_year}"

    # Retrieve the name for the specified period.
    if query_period in archimandrites_data:
        name = archimandrites_data[query_period]
        print(f"The archimandrite of the Pskov-Caves Monastery from {start_year} to {end_year} was: {name}")
    else:
        print(f"No data found for the period {query_period}.")

find_archimandrite()