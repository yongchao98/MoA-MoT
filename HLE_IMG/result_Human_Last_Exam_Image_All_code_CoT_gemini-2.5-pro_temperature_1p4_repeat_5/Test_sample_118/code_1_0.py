def solve_zurich_path():
    """
    This function calculates and prints the best path from Billoweg to Laubiweg
    based on the provided Zurich public transport map and criteria.
    The path is determined by manual analysis of the map.
    """
    # Details of the optimal path found
    start_symbol = "S"
    end_symbol = "E"
    start_line = 72
    exchange_station = "Hertersteg"
    end_line = 32
    total_stations = 8

    # The final path description string
    path_description = (
        f"{start_symbol} - {start_line} - {exchange_station} - "
        f"{end_line} - {end_symbol}; {total_stations}"
    )

    print(path_description)

solve_zurich_path()