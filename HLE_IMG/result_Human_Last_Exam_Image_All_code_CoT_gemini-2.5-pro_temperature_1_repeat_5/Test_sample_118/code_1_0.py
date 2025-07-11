def find_best_path():
    """
    Finds the best path from Billoweg to Laubiweg based on given criteria.
    """
    # Step 1: Model the relevant parts of the network
    # Data is transcribed from the provided map for the relevant lines.
    network = {
        '7': [
            'Billoweg', 'Bhf. Wollishofen/Staubstr.', 'Renggerstr.', 'Morgental', 
            'Poststrasse', 'Brunaustr.', 'Museum Rietberg', 'Bahnhof Enge', 
            'Tunnelstr.', 'Stockerstr.', 'Paradeplatz', 'Rennweg', 
            'Bahnhofstrasse/HB', 'Bahnhofquai/HB', 'Central', 'Haldenegg', 
            'Sonneggstr.', 'Ottikerstr.', 'Schaffhauserplatz', 'Guggachstr.', 
            'Milchbuck'
        ],
        '66': [
            'Bahnhof Enge', 'Brunaustr.', 'Hügelstr.', 'Sihlcity Nord', 
            'Saalsporthalle/Sihlcity', 'Laubiweg', 'Brunau'
        ],
        '72': [
            'Milchbuck', 'Guggachstr.', 'Berninaplatz', 'Brunnenhof', 
            'Bucheggplatz', 'Lehen-Steig', 'Escher-gutweg', 'Rebbergsteig', 
            'Wipkingerplatz', 'Rosengartenstr.', 'Hardplatz', 'Albisriederplatz', 
            'Zypressenstr.', 'Lochergut', 'Kalkbreite/Bhf. Wiedikon', 
            'Schmiede Wiedikon', 'Manessestrasse', 'Giesshübel', 'Sihlcity Nord', 
            'Saalsporthalle/Sihlcity', 'Laubiweg', 'Brunau', 'Morgental'
        ]
    }
    
    start_station = 'Billoweg'
    end_station = 'Laubiweg'

    # Helper function to get the number of stations in a segment
    def get_segment_station_count(line, station1, station2):
        try:
            stations_on_line = network[line]
            idx1 = stations_on_line.index(station1)
            idx2 = stations_on_line.index(station2)
            return abs(idx1 - idx2) + 1
        except (KeyError, ValueError):
            return float('inf')

    # Step 2: Search for paths, starting with the lowest number of exchanges.
    # We will search for 0-exchange, then 1-exchange paths.

    all_paths = []
    
    # Identify lines serving start and end stations
    start_lines = [line for line, stations in network.items() if start_station in stations]
    end_lines = [line for line, stations in network.items() if end_station in stations]

    # --- Search for 0-exchange paths ---
    # This would be a direct connection on a single line.
    common_lines = set(start_lines).intersection(set(end_lines))
    if common_lines:
        for line in common_lines:
            station_count = get_segment_station_count(line, start_station, end_station)
            path_details = {
                'exchanges': 0,
                'stations': station_count,
                'path_str': f"S - {line} - E",
                'exchange_pos': station_count # For sorting, exchange is at the end
            }
            all_paths.append(path_details)
    
    # If we found a direct path, we don't need to look for paths with more exchanges.
    # In this problem, no 0-exchange path exists, so we proceed.
    
    # --- Search for 1-exchange paths ---
    if not all_paths:
        for start_line in start_lines:
            for end_line in end_lines:
                if start_line == end_line:
                    continue
                
                # Find all possible exchange stations
                start_line_stations = set(network[start_line])
                end_line_stations = set(network[end_line])
                exchange_stations = start_line_stations.intersection(end_line_stations)

                for exchange in exchange_stations:
                    # Calculate path properties
                    len1 = get_segment_station_count(start_line, start_station, exchange)
                    len2 = get_segment_station_count(end_line, exchange, end_station)
                    
                    if len1 != float('inf') and len2 != float('inf'):
                        total_stations = len1 + len2 - 1
                        path_details = {
                            'exchanges': 1,
                            'stations': total_stations,
                            'path_str': f"S - {start_line} - {exchange} - {end_line} - E",
                            'exchange_pos': len1 # Criterion C: latest possible exchange
                        }
                        all_paths.append(path_details)

    # Step 3: Select the best path based on the criteria
    # Sort by: 1. exchanges (asc), 2. stations (asc), 3. exchange_pos (desc)
    if not all_paths:
        print("No path found.")
        return

    all_paths.sort(key=lambda p: (p['exchanges'], p['stations'], -p['exchange_pos']))
    
    best_path = all_paths[0]

    # Step 4: Format and print the final answer
    # The prompt asks to output each number in the final equation.
    # The final string will contain all required numbers.
    final_answer_string = f"{best_path['path_str']}; {best_path['stations']}"
    print(final_answer_string)


if __name__ == '__main__':
    find_best_path()
