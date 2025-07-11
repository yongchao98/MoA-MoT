def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided Zurich public transport map and a specific set of criteria.

    The criteria in order of importance are:
    1. Minimum number of exchanges.
    2. Minimum number of stations.
    3. Latest possible exchanges.

    The analysis identified that the best path requires 2 exchanges and has 16 stations.
    The chosen path is preferred because its first exchange occurs at the 10th station,
    which is later than the alternative path with 16 stations (which exchanges at the 3rd station).
    """

    # The optimal path found through analysis:
    # Start (S) at Billoweg.
    line1 = 72
    exchange1 = "Schmiede Wiedikon"
    # Take Bus 72 to Schmiede Wiedikon.
    
    line2 = 9
    exchange2 = "Sihlpost/HB"
    # Change to Tram 9, travel to Sihlpost/HB.
    
    line3 = 15
    # Change to Tram 15, travel to Laubiweg (E).
    
    station_count = 16
    # The total number of unique stations on the path is 16.

    # Construct the final output string in the specified format.
    # The string includes the line numbers (72, 9, 15) and the total station count (16).
    final_path = f"S - {line1} - {exchange1} - {line2} - {exchange2} - {line3} - E; {station_count}"

    print(final_path)

solve_zurich_path()