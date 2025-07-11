def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided criteria.
    """
    start_station = "S"
    end_station = "E"
    
    # The optimal path found after analyzing the map
    line1 = 7
    exchange_station = "Bahnhofquai/HB"
    line2 = 15
    total_stations = 19
    
    # Constructing the output string in the specified format
    path_description = (
        f"{start_station} - {line1} - {exchange_station} - "
        f"{line2} - {end_station}; {total_stations}"
    )
    
    print(path_description)

solve_zurich_path()