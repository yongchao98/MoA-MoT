def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided criteria and map of Zurich's public transport system.
    """
    start_station = "Billoweg"
    end_station = "Laubiweg"

    # The optimal path determined by applying the given criteria:
    # 1. Start at Billoweg on Line 7.
    # 2. Travel 18 stations to the exchange point.
    first_line = 7
    exchange_station = "Schaffhauserplatz"

    # 3. Switch to Line 33.
    # 4. Travel 2 stations to the destination, Laubiweg.
    second_line = 33

    # 5. Calculate the total number of stations.
    # (Stations on Leg 1) + (Stations on Leg 2) - 1 (for the exchange station counted twice)
    # Leg 1: Billoweg -> Schaffhauserplatz on Line 7 is 18 stations.
    # Leg 2: Schaffhauserplatz -> Laubiweg on Line 33 is 2 stations.
    total_stations = 18 + 2 - 1
    
    # Format the output string as requested.
    # S - 7 - Schaffhauserplatz - 33 - E; 19
    path_description = f"S - {first_line} - {exchange_station} - {second_line} - E; {total_stations}"
    
    print(path_description)

solve_zurich_path()