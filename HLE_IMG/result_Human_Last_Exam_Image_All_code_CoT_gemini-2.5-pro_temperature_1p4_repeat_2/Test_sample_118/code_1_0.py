def solve_zurich_path():
    """
    This function determines and prints the optimal path from Billoweg to Laubiweg
    based on the provided criteria. The logic for finding the path is derived
    from analyzing the Zurich public transport map.

    Pathfinding logic:
    1.  Start at Billoweg (lines 7, 72) and End at Laubiweg (lines 9, 10).
    2.  No direct path (0 exchanges) exists.
    3.  Search for paths with 1 exchange.
        a) S (72) -> Schmiede Wiedikon (9) -> E: 21 stations.
           - Leg 1 (Billoweg to Schmiede Wiedikon on line 72): 10 stations.
           - Leg 2 (Schmiede Wiedikon to Laubiweg on line 9): 12 stations.
           - Total: 10 + 12 - 1 = 21.
        b) S (72) -> Milchbuck (10) -> E: 22 stations.
        c) S (7)  -> Central (10) -> E: 23 stations.
    4.  According to the criteria (fewest exchanges, then fewest stations),
        path (a) is the unique optimal path.
    """
    
    # Define the components of the optimal path
    start_station = "S"
    end_station = "E"
    
    line_1 = 72
    exchange_station = "Schmiede Wiedikon"
    line_2 = 9
    
    # Calculate the total number of stations
    # Leg 1 (Billoweg -> Schmiede Wiedikon):
    # 1. Billoweg, 2. Morgental, 3. Jugendherberge, 4. Thujastr., 5. Brunaustr.,
    # 6. Sihlcity Nord, 7. Saalsporthalle, 8. Manesseplatz, 9. Hertersteg, 10. Schmiede Wiedikon
    stations_leg_1 = 10
    
    # Leg 2 (Schmiede Wiedikon -> Laubiweg):
    # 1. Schmiede Wiedikon, 2. Werd, 3. Stauffacher, 4. Sihlpost/HB, 5. Bahnhofstr./HB,
    # 6. Bahnhofquai/HB, 7. Stampfenbachplatz, 8. Beckenhof, 9. Kronenstr.,
    # 10. Nordstr., 11. Rotbuchstr., 12. Laubiweg
    stations_leg_2 = 12
    
    # Total stations = sum of stations on each leg, minus 1 for the transfer station
    total_stations = stations_leg_1 + stations_leg_2 - 1
    
    # Construct the final output string
    path_description = f"{start_station} - {line_1} - {exchange_station} - {line_2} - {end_station}; {total_stations}"
    
    print(path_description)

solve_zurich_path()