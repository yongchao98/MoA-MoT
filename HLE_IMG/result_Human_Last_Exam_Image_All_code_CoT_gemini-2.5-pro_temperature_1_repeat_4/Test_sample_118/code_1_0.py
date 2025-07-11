def solve_zurich_path():
    """
    Finds the optimal path from Billoweg to Laubiweg based on the given criteria.
    """
    # Step 1: Model the network with manually extracted line data.
    # The lists represent the sequence of stations on each line.
    lines = {
        "72": [
            "Morgental", "Billoweg", "Jugendherberge", "Thujastr", "Brunaustr", "Hügelstr",
            "Saalsporthalle", "Bhf. Enge/Bederstr", "Tunnelstr", "Stockerstr", "Paradeplatz",
            "Helmhaus", "Bellevue", "Opernhaus", "Kreuzplatz", "Signaustr", "Englischviertelstr",
            "Sprecherstr", "Bhf. Stadelhofen", "Kunsthaus", "Neumarkt", "Central", "Haldenegg",
            "Sonneggstr", "Scheuchzerstr", "Ottikerstr", "Röslistr", "Schaffhauserplatz",
            "Guggachstr", "Milchbuck"
        ],
        "9": [
            "Triemli", "Heuried", "Goldbrunnenplatz", "Schmiede Wiedikon", "Kalkbreite/Bhf Wiedikon",
            "Stauffacher", "Sihlstr", "Kantonalbank", "Paradeplatz", "Rennweg", "Bahnhofquai/HB",
            "Stampfenbachplatz", "Kronenstr", "Schaffhauserplatz", "Guggachstr", "Laubiweg", "Milchbuck"
        ],
        "15": [
            "Klusplatz", "Kreuzplatz", "Bhf. Stadelhofen", "Bellevue", "Paradeplatz", "Löwenplatz",
            "Bahnhofplatz/HB", "Central", "Haldenegg", "Sonneggstr", "Letzistr", "Röslistr",
            "Ottikerstr", "Schaffhauserplatz", "Guggachstr", "Laubiweg", "Bucheggplatz"
        ],
        "32": [
            "Strassenverkehrsamt", "Albisgütli", "Siemens", "Hubertus", "Langgrütstr", "Im Gut",
            "Heuried", "Goldbrunnenplatz", "Schmiede Wiedikon", "Kalkbreite/Bhf Wiedikon",
            "Helvetiaplatz", "Militär-/Langstr", "Quellenstr", "Limmatplatz", "Rotbuchstr",
            "Radiostudio", "Bucheggplatz", "Laubiweg", "Guggachstr", "Waidbadstr", "Wartau",
            "Winzerhalde", "Meierhofplatz", "Schützenhaus Höngg", "Friedhof Hönggerberg",
            "Zehntenhausplatz", "Chrüzächer", "Holzerhurd"
        ]
    }

    start_station = "Billoweg"
    end_station = "Laubiweg"

    # Step 2: Identify potential paths (focus on 1-exchange paths as per Criterion A).
    start_lines = [name for name, stations in lines.items() if start_station in stations]
    end_lines = [name for name, stations in lines.items() if end_station in stations]
    
    possible_paths = []

    for start_line_name in start_lines:
        for end_line_name in end_lines:
            if start_line_name == end_line_name:
                continue

            start_line_stations = lines[start_line_name]
            end_line_stations = lines[end_line_name]

            # Find all intersection points between the two lines.
            intersections = set(start_line_stations) & set(end_line_stations)

            for exchange_station in intersections:
                # Step 3: Evaluate paths by station count (Criterion B).
                try:
                    # Calculate stations on the first leg of the journey.
                    start_idx = start_line_stations.index(start_station)
                    exchange_idx_1 = start_line_stations.index(exchange_station)
                    stations_leg1 = abs(start_idx - exchange_idx_1) + 1

                    # Calculate stations on the second leg.
                    exchange_idx_2 = end_line_stations.index(exchange_station)
                    end_idx = end_line_stations.index(end_station)
                    stations_leg2 = abs(exchange_idx_2 - end_idx) + 1

                    total_stations = stations_leg1 + stations_leg2 - 1
                    
                    path_details = {
                        "line1": start_line_name,
                        "exchange": exchange_station,
                        "line2": end_line_name,
                        "stations": total_stations
                    }
                    possible_paths.append(path_details)
                except ValueError:
                    # This station does not exist on this line, skip.
                    continue

    # Step 4: Select the optimal path.
    # Sort paths by the number of stations (fewer is better).
    if not possible_paths:
        print("No path with one exchange found.")
        return

    # Criterion A (fewer exchanges) is met by only considering 1-exchange paths.
    # Criterion B (fewer stations) is used for sorting.
    possible_paths.sort(key=lambda p: p['stations'])
    best_path = possible_paths[0]

    # Step 5: Format and print the output.
    line1 = best_path['line1']
    exchange_station = best_path['exchange']
    line2 = best_path['line2']
    station_count = best_path['stations']

    # The final output string includes all the required numbers.
    result_string = f"S - {line1} - {exchange_station} - {line2} - E; {station_count}"
    print(result_string)

solve_zurich_path()
<<<S - 72 - Paradeplatz - 9 - E; 17>>>