def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided criteria and map.
    """
    
    # Path details derived from map analysis
    start_symbol = "S"
    first_line = 72
    exchange_station = "Milchbuck"
    second_line = 32
    end_symbol = "E"
    
    # Station count calculation
    # Segment 1 (Line 72: Billoweg -> Milchbuck): 19 stations
    # Segment 2 (Line 32: Milchbuck -> Laubiweg): 17 stations
    # Total = 19 + 17 - 1 (exchange station counted once)
    total_stations = 35
    
    # Format and print the final answer
    # The final output string includes all the required numbers.
    final_answer = f"{start_symbol} - {first_line} - {exchange_station} - {second_line} - {end_symbol}; {total_stations}"
    print(final_answer)

solve_zurich_path()
<<<S - 72 - Milchbuck - 32 - E; 35>>>