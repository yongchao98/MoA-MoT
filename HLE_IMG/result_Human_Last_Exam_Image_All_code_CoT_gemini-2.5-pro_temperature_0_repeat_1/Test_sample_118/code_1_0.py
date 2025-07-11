def find_best_path():
    """
    This function analyzes potential paths from Billoweg to Laubiweg based on
    pre-determined data from the transport map and finds the optimal one
    according to the specified criteria.
    """
    # Data representing the possible paths with 1 exchange.
    # 'stations_leg1' is the number of stations before the exchange.
    # 'total_stations' is the total number of stations on the path.
    paths = [
        {
            "line1": 72, "exchange1": "Schmiede Wiedikon", "line2": 9,
            "exchanges": 1, "stations_leg1": 13, "total_stations": 23
        },
        {
            "line1": 72, "exchange1": "Milchbuck", "line2": 9,
            "exchanges": 1, "stations_leg1": 21, "total_stations": 26
        },
        {
            "line1": 72, "exchange1": "Bucheggplatz", "line2": 15,
            "exchanges": 1, "stations_leg1": 20, "total_stations": 22
        },
        {
            "line1": 72, "exchange1": "Schmiede Wiedikon", "line2": 32,
            "exchanges": 1, "stations_leg1": 13, "total_stations": 21
        },
        {
            "line1": 72, "exchange1": "Bucheggplatz", "line2": 32,
            "exchanges": 1, "stations_leg1": 20, "total_stations": 22
        }
    ]

    # Sort the paths based on the criteria:
    # 1. Number of exchanges (ascending).
    # 2. Total number of stations (ascending).
    # 3. Number of stations in the first leg (descending, for latest exchange).
    # The negative sign on 'stations_leg1' achieves descending order for the third criterion.
    sorted_paths = sorted(paths, key=lambda p: (p['exchanges'], p['total_stations'], -p['stations_leg1']))

    # The best path is the first one in the sorted list.
    best_path = sorted_paths[0]

    # Format the output string as required.
    output_string = (f"S - {best_path['line1']} - {best_path['exchange1']} - "
                     f"{best_path['line2']} - E; {best_path['total_stations']}")

    print(output_string)

find_best_path()