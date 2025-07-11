def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided criteria.
    """
    # 1. Define the components of the optimal path based on map analysis.
    # The optimal path has 1 exchange, which is better than 2 or more.
    # Comparing 1-exchange paths, the one via Guggachstr has 15 stations and a late exchange,
    # making it optimal according to the given criteria.

    start_symbol = "S"
    end_symbol = "E"

    # First leg of the journey
    line1 = 72

    # Exchange station
    exchange_station = "Guggachstr"

    # Second leg of the journey
    line2 = 15

    # 2. Calculate the total number of stations in the optimal path.
    # The stations are counted inclusively from start to end.
    # S(Billoweg) -> Morgental -> Jugendherberge -> Brunaustr. -> Bahnhof Enge -> Tunnelstr. -> Stockerstr.
    # -> Sihlstr./Selnau -> Löwenplatz -> Bahnhofplatz/HB -> Stampfenbachplatz -> Beckenhof
    # -> Schaffhauserplatz -> Guggachstr (Exchange) -> E(Laubiweg)
    total_stations = 15

    # 3. Format the final output string as per the specified format.
    # Format: S - <line> - <station> - <line> - E; <count>
    # Note: the prompt asks to output each number in the final equation.
    path_description = f"{start_symbol} - {line1} - {exchange_station} - {line2} - {end_symbol}; {total_stations}"

    print(path_description)

solve_zurich_path()
<<<S - 72 - Guggachstr - 15 - E; 15>>>