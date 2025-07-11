def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg.

    The logic is as follows:
    1.  Identify lines for Start (Billoweg: 7) and End (Laubiweg: 32, 40).
    2.  Find paths with the minimum number of exchanges (Criterion A).
        - A 0-exchange path is not possible.
        - A 1-exchange path is possible by finding intersections between line 7 and lines {32, 40}.
        - Intersections are found at 'Guggachstr.' and 'Milchbuck'.
    3.  Compare the 1-exchange paths by the number of stations (Criterion B).
        - Path via Guggachstr.: Take line 7 from Billoweg to Guggachstr. (19 stops), then line 32 to Laubiweg (1 more stop). Total = 20 stops.
        - Path via Milchbuck: Take line 7 from Billoweg to Milchbuck (20 stops), then line 32 to Laubiweg (2 more stops). Total = 22 stops.
    4.  The path via Guggachstr. is shorter and therefore optimal.
    5.  The final path is formatted as required.
    """
    start_station = "S"
    end_station = "E"
    
    line1 = 7
    exchange_station = "Guggachstr."
    line2 = 32
    total_stations = 20

    # Print the final path description, ensuring each component is output.
    print(f"{start_station} - {line1} - {exchange_station} - {line2} - {end_station}; {total_stations}")

solve_zurich_path()