def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided map and criteria.
    """
    # Path details deduced from the map analysis
    start_station = "S"
    end_station = "E"
    first_line = 7
    exchange_station = "Guggachstr."
    second_line = 15
    total_stations = 19

    # Formatting the output string as per the user's request
    # The output format is: S - <line1> - <exchange1> - <line2> - E; <station_count>
    path_description = (
        f"{start_station} - {first_line} - {exchange_station} - "
        f"{second_line} - {end_station}; {total_stations}"
    )

    print(path_description)

solve_zurich_path()