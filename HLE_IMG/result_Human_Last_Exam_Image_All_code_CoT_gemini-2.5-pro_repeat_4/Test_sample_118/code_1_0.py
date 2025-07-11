def solve_zurich_tram_path():
    """
    This function determines the optimal path from Billoweg to Laubiweg
    based on the provided criteria and prints the result.
    The logic is derived from analyzing the transport map as detailed in the thought process.
    """
    # Based on the step-by-step analysis of the map:
    # 1. The path with the minimum number of exchanges (1) is chosen.
    # 2. Between the 1-exchange options, the one with the minimum number of stations is chosen.
    #    - Path via Schmiede Wiedikon: 12 stations
    #    - Path via Escher-Wyss-Platz: 21 stations
    #    The path via Schmiede Wiedikon is optimal.
    # 3. The ambiguity at Schmiede Wiedikon (take Tram 9 or 14) is resolved by choosing the lower line number to ensure uniqueness.
    
    start_station = "S"
    end_station = "E"
    
    line1 = 72
    exchange_station1 = "Schmiede Wiedikon"
    line2 = 9
    total_stations = 12
    
    # Format the output as specified in the problem description.
    # S - [line1] - [exchange1] - [line2] - E; [station_count]
    final_path_description = f"{start_station} - {line1} - {exchange_station1} - {line2} - {end_station}; {total_stations}"
    
    print(final_path_description)

solve_zurich_tram_path()