def solve_zurich_path():
    """
    This function provides the solution to the Zurich tram pathfinding problem
    from Billoweg to Laubiweg based on a step-by-step analysis of the provided map
    and criteria.

    The reasoning is as follows:
    1.  **Criterion A (Fewer Exchanges):** A path with 1 exchange is optimal, as no 0-exchange path exists.
    2.  **Criterion B (Fewer Stations):** We compare the two possible 1-exchange paths:
        *   **Path 1 (via Hardbrücke):** Take Bus 72 to Hardbrücke, then switch to Bus 32. This path involves 22 stations.
        *   **Path 2 (via Bucheggplatz):** Take Bus 72 to Bucheggplatz, then switch to another line. This path involves 20 stations.
        Path 2 is superior as it has fewer stations.
    3.  **Uniqueness and Final Choice:** The optimal path switches at Bucheggplatz. The final segment from Bucheggplatz to Laubiweg is served by three lines: Tram 11, Tram 15, and Bus 32. All three choices result in a path with the same number of exchanges (1), the same number of stations (20), and the same exchange point. Since the problem states a unique path exists, a tie-breaking rule is needed. A standard logical tie-breaker is to select the line with the lowest number. Out of {11, 15, 32}, the lowest is 11.
    4.  **Final Path Description:**
        *   Start at Billoweg (S).
        *   Take line 72 to Bucheggplatz.
        *   Switch to line 11.
        *   Arrive at Laubiweg (E).
        *   Total stations: 20.
    """
    start_symbol = "S"
    first_line_number = 72
    exchange_station_name = "Bucheggplatz"
    second_line_number = 11
    end_symbol = "E"
    total_stations = 20

    # The final answer is formatted as requested.
    # The instruction "you still need to output each number in the final equation!" is interpreted
    # as ensuring the numbers are present and clear in the final printed string.
    print(f"{start_symbol} - {first_line_number} - {exchange_station_name} - {second_line_number} - {end_symbol}; {total_stations}")

solve_zurich_path()
<<<S - 72 - Bucheggplatz - 11 - E; 20>>>