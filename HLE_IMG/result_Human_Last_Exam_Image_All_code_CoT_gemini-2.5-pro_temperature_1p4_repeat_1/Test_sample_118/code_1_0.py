def solve_zurich_path():
    """
    This function formats and prints the optimal path from Billoweg to Laubiweg
    based on the analysis of the provided transport map and criteria.
    """
    
    # Path details determined by map analysis
    start_point = "S"
    end_point = "E"
    
    # First leg of the journey
    line1 = 72
    exchange_station = "Bucheggplatz"
    
    # Second leg of the journey
    line2 = 15
    
    # Total stations travelled
    total_stations = 18
    
    # Constructing the final path string according to the specified format
    # S - 72 - Bucheggplatz - 15 - E; 18
    final_path = (
        f"{start_point} - {line1} - {exchange_station} - "
        f"{line2} - {end_point}; {total_stations}"
    )
    
    print(final_path)

solve_zurich_path()