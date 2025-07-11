import collections

def find_optimal_path():
    """
    This function determines the best path from Billoweg to Laubiweg by analyzing
    and comparing potential routes based on a set of hierarchical criteria.
    """
    # The analysis of the Zurich transport map reveals that the optimal path must have
    # one exchange, transferring from Line 72 to Line 17. There are four possible
    # stations for this exchange.
    
    # Manually transcribed station data for the relevant routes
    line_72_stations_from_start = [
        "Billoweg", "Brunaustr.", "Hügelstr.", "Sihlcity Nord", "Giesshübel",
        "Waffenplatzstr.", "Bahnhof Selnau", "Tunnelstr.", "Stockerstr.",
        "Paradeplatz", "Rennweg", "Bahnhofstr./HB", "Bahnhofquai/HB", "Central"
    ]

    line_17_stations_to_end = [
        "Central", "Bahnhofquai/HB", "Bahnhofstr./HB", "Rennweg", "Sihlpost",
        "Militär-/Langstr.", "Albisriederplatz", "Hubertus", "Im Gut",
        "Heuried", "Laubegg", "Laubiweg"
    ]

    # The four stations where a transfer from Line 72 to Line 17 is possible
    exchange_points = ["Rennweg", "Bahnhofstr./HB", "Bahnhofquai/HB", "Central"]

    candidate_paths = []
    Path = collections.namedtuple('Path', ['line_1', 'exchange_station', 'line_2', 'num_exchanges', 'num_stations', 'lateness'])

    for station in exchange_points:
        # Determine the stations on the first leg of the journey
        leg1_end_index = line_72_stations_from_start.index(station)
        leg1_stations = line_72_stations_from_start[:leg1_end_index + 1]
        
        # Determine the stations on the second leg
        leg2_start_index = line_17_stations_to_end.index(station)
        leg2_stations = line_17_stations_to_end[leg2_start_index:]
        
        # Criterion B: Calculate the total number of stations
        total_stations = len(leg1_stations) + len(leg2_stations) - 1
        
        # Criterion C: Measure the "lateness" of the exchange
        lateness_score = len(leg1_stations)

        candidate_paths.append(Path(
            line_1=72,
            exchange_station=station,
            line_2=17,
            num_exchanges=1,
            num_stations=total_stations,
            lateness=lateness_score
        ))

    # Sort paths by criteria: A) exchanges (asc), B) stations (asc), C) lateness (desc)
    sorted_paths = sorted(candidate_paths, key=lambda p: (p.num_exchanges, p.num_stations, -p.lateness))
    
    best_path = sorted_paths[0]

    # Print the description of the best path in the specified format
    print(f"S - {best_path.line_1} - {best_path.exchange_station} - {best_path.line_2} - E; {best_path.num_stations}")

find_optimal_path()