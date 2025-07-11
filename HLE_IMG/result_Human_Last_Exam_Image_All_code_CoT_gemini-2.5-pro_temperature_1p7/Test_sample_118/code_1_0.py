def find_optimal_path():
    """
    This function describes the optimal path from Billoweg to Laubiweg based on a manual analysis
    of the Zurich tram map and the provided criteria.

    Analysis Summary:
    1.  Criterion A (Exchanges): No 0-exchange path exists. We search for 1-exchange paths.
    2.  Possible 1-Exchange Paths:
        a) From Billoweg via line 72, exchange at "Schmiede Wiedikon" to line 32, then to Laubiweg.
        b) From Billoweg via line 72, exchange at "Bucheggplatz" to lines 15, 32, or 40, then to Laubiweg.
    3.  Criterion B (Stations):
        - Path (a) via Schmiede Wiedikon:
            - Leg 1 (Billoweg -> Schmiede Wiedikon on line 72): 8 stations.
            - Leg 2 (Schmiede Wiedikon -> Laubiweg on line 32): 11 stations.
            - Total: 8 + 11 - 1 = 18 stations.
        - Path (b) via Bucheggplatz:
            - Leg 1 (Billoweg -> Bucheggplatz on line 72): 19 stations.
            - Leg 2 (Bucheggplatz -> Laubiweg on lines 15/32/40): 2 stations.
            - Total: 19 + 2 - 1 = 20 stations.
    4.  Conclusion: The path via "Schmiede Wiedikon" is optimal because it has the same number of exchanges (1) but fewer total stations (18 vs. 20).
    """

    # Path Details
    start_station_symbol = "S"
    first_line = 72
    exchange_station = "Schmiede Wiedikon"
    second_line = 32
    end_station_symbol = "E"
    total_stations = 18

    # Format the final path description as per the instructions
    # S - 72 - Schmiede Wiedikon - 32 - E; 18
    final_path = f"{start_station_symbol} - {first_line} - {exchange_station} - {second_line} - {end_station_symbol}; {total_stations}"
    
    print(final_path)

find_optimal_path()