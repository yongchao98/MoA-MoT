def solve_zurich_path():
    """
    This function determines the optimal path from Billoweg to Laubiweg
    based on the provided map and criteria, then prints the result.
    """
    
    # Path details based on analysis of the Zurich transport map
    start_station = "S"
    end_station = "E"
    
    first_line = 72
    exchange_station = "Saalsporthalle"
    second_line = 5
    total_stations = 7
    
    # Construct the final path description string
    path_description = (
        f"{start_station} - {first_line} - {exchange_station} - "
        f"{second_line} - {end_station}; {total_stations}"
    )
    
    print(path_description)

solve_zurich_path()