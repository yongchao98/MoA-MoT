def solve_zurich_path():
    """
    This function determines the optimal path from Billoweg to Laubiweg
    based on the provided map and criteria, and prints the result in the specified format.
    """
    start_station = "S"
    end_station = "E"
    
    # Based on the step-by-step analysis of the map and criteria.
    # Criterion A: Minimum exchanges (1) is preferred.
    # Criterion B: Minimum stations for a 1-exchange path is 21.
    # Criterion C: The path with the latest exchange point (Bahnhofplatz/HB) is preferred.
    # Tie-breaker (9 vs 14): The station info on the map at Bahnhofplatz/HB shows Line 14, not 9.
    
    initial_line = 7
    exchange_station = "Bahnhofplatz/HB"
    final_line = 14
    total_stations = 21 # 14 stops on Line 7 + 8 stops on Line 14 - 1 for the transfer station.

    path_description = f"{start_station} - {initial_line} - {exchange_station} - {final_line} - {end_station}; {total_stations}"
    
    print(path_description)

solve_zurich_path()