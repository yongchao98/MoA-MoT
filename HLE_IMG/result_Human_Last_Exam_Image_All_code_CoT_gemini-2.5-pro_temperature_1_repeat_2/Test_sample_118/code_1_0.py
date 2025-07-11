def find_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on a pre-analysis of the Zurich public transport map and the given criteria.

    The criteria in descending order of importance are:
    A) Fewer exchanges
    B) Fewer stations
    C) Latest possible exchanges for paths that are otherwise identical

    Analysis Summary:
    1.  A 0-exchange path is not possible. The minimum is 1 exchange.
    2.  Two 1-exchange paths have the minimum number of stations (18):
        - Path A: S --(72)--> Bahnhofquai/HB --(15)--> E; 18 stations
        - Path B: S --(72)--> Stampfenbachplatz --(15)--> E; 18 stations
    3.  Comparing these two paths with Criterion C, the exchange at Stampfenbachplatz
        occurs later in the journey on line 72 than the one at Bahnhofquai/HB.
        Therefore, Path B is the unique optimal path.
    """

    # Details of the optimal path
    start_symbol = "S"
    end_symbol = "E"
    first_line = 72
    exchange_station = "Stampfenbachplatz"
    second_line = 15
    total_stations = 18

    # Construct the final output string as per the required format
    # The string includes the line numbers (72, 15) and the total station count (18).
    result = f"{start_symbol} - {first_line} - {exchange_station} - {second_line} - {end_symbol}; {total_stations}"

    print(result)

find_zurich_path()