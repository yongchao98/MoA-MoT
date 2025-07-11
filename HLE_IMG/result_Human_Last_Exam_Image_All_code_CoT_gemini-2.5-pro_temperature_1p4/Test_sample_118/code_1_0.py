def solve_zurich_path():
    """
    This script determines the optimal path from Billoweg to Laubiweg
    based on the provided criteria and formats the result.
    """
    
    # Path analysis based on the map and criteria.
    # Criterion A: Minimize exchanges.
    # A direct path on line 72 exists, which means 0 exchanges. This is optimal.
    
    # Details of the optimal path
    start_symbol = "S"
    end_symbol = "E"
    line_number = 72
    
    # Criterion B: Minimize stations.
    # The stations on the path are:
    # 1. Billoweg (Start)
    # 2. Brunaustr.
    # 3. Jugendherberge
    # 4. Besenrainstr.
    # 5. Sihlcity
    # 6. Saalsporthalle
    # 7. Laubiweg (End)
    station_count = 7
    
    # Format the final output string as per the problem description.
    # The format is "S - 72 - E; 7"
    # This shows the start, the line taken, the end, and the total station count.
    
    path_description = f"{start_symbol} - {line_number} - {end_symbol}"
    final_answer = f"{path_description}; {station_count}"
    
    print(final_answer)

solve_zurich_path()