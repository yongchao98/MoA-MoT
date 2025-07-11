def solve_zurich_path():
    """
    This function calculates and prints the optimal path from Billoweg to Laubiweg
    based on the provided Zurich public transport map and criteria.
    """
    
    # Path details derived from map analysis
    start_station = "S"
    end_station = "E"
    
    # Path segment 1
    line1 = 72
    
    # Exchange
    exchange_station_1 = "Hubertus"
    
    # Path segment 2
    line2 = 89
    
    # Total number of stations traversed
    total_stations = 11
    
    # Construct the final path description string
    path_description = f"{start_station} - {line1} - {exchange_station_1} - {line2} - {end_station}; {total_stations}"
    
    print(path_description)
    
    # The final answer in the required format
    final_answer = f"<<<{path_description}>>>"
    print(final_answer)

solve_zurich_path()