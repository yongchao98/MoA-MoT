def solve_zurich_path():
    """
    This function prints the solution for the Zurich tram path problem.
    The path is determined by following the specified criteria:
    1. Minimize exchanges.
    2. Minimize stations.
    3. Exchange as late as possible.
    
    The optimal path from Billoweg to Laubiweg is:
    - Start at Billoweg, take bus 72 to Wipkingerplatz. (14 stops)
    - At Wipkingerplatz, switch to tram 9.
    - Take tram 9 to Laubiweg. (4 stops)
    - Total stops = 14 + 4 - 1 = 17.
    - An alternative path via Sihlpost/HB also has 17 stops but the exchange is earlier (after 7 stops), so it is not optimal.
    """
    start = "S"
    end = "E"
    
    first_line = 72
    exchange_station = "Wipkingerplatz"
    second_line = 9
    
    total_stations = 17
    
    # Printing each component of the final answer string as requested
    print(f"{start} - {first_line} - {exchange_station} - {second_line} - {end}; {total_stations}")

solve_zurich_path()