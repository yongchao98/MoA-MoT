def solve_zurich_path():
    """
    This function determines the best path from Billoweg to Laubiweg
    based on the provided criteria and prints the result in the specified format.
    """
    # Path analysis details:
    # 1. Start: Billoweg (Line 72)
    # 2. End: Laubiweg (Lines 9, 14, 32)
    # 3. Minimum exchanges: 1. A direct path on a single line is not possible.
    #    The best exchange point is 'Schmiede Wiedikon', where Bus 72 intersects Trams 9 and 14.
    #    Other exchange points result in a much higher station count.
    # 4. Minimum stations: The path via Schmiede Wiedikon has 7 stations.
    #    Path: Billoweg -> Brunaustr. -> Giesshübel (Sihlcity Nord) -> Waffenplatzstr. ->
    #    Schmiede Wiedikon -> Goldbrunnenplatz -> Laubiweg.
    # 5. Tie-breaking: At Schmiede Wiedikon, we can switch to either Tram 9 or Tram 14.
    #    Both paths are identical in terms of exchanges, station count, and exchange location.
    #    To find the unique path promised, a tie-breaker is assumed: choose the lower line number.
    #    Since 9 < 14, we choose Tram line 9.

    start_station = "S"
    end_station = "E"
    
    first_line = 72
    exchange_station = "Schmiede Wiedikon"
    second_line = 9
    total_stations = 7

    # Constructing the final path description string
    # The final equation requires all the numbers to be present.
    final_path = f"{start_station} - {first_line} - {exchange_station} - {second_line} - {end_station}; {total_stations}"
    
    print(final_path)

solve_zurich_path()