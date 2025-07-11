def count_mlb_triple_crowns():
    """
    This function calculates and displays the number of times the MLB batting Triple Crown has been won.
    """
    # Data source: A comprehensive list of MLB Batting Triple Crown winners.
    # A player wins the Triple Crown by leading their league in batting average, home runs, and RBIs.
    triple_crown_instances = [
        # American League Winners
        {"year": 1909, "name": "Ty Cobb", "league": "AL"},
        {"year": 1933, "name": "Jimmie Foxx", "league": "AL"},
        {"year": 1934, "name": "Lou Gehrig", "league": "AL"},
        {"year": 1942, "name": "Ted Williams", "league": "AL"},
        {"year": 1947, "name": "Ted Williams", "league": "AL"},
        {"year": 1956, "name": "Mickey Mantle", "league": "AL"},
        {"year": 1966, "name": "Frank Robinson", "league": "AL"},
        {"year": 1967, "name": "Carl Yastrzemski", "league": "AL"},
        {"year": 2012, "name": "Miguel Cabrera", "league": "AL"},
        
        # National League Winners
        {"year": 1878, "name": "Paul Hines", "league": "NL"},
        {"year": 1894, "name": "Hugh Duffy", "league": "NL"},
        {"year": 1922, "name": "Rogers Hornsby", "league": "NL"},
        {"year": 1925, "name": "Rogers Hornsby", "league": "NL"},
        {"year": 1933, "name": "Chuck Klein", "league": "NL"},
        {"year": 1937, "name": "Joe Medwick", "league": "NL"},
    ]

    al_count = sum(1 for winner in triple_crown_instances if winner['league'] == 'AL')
    nl_count = sum(1 for winner in triple_crown_instances if winner['league'] == 'NL')
    total_count = len(triple_crown_instances)

    print("Counting the MLB Batting Triple Crown winners by league:")
    print(f"American League winners: {al_count}")
    print(f"National League winners: {nl_count}")
    print("\nTotal number of Triple Crown wins in MLB history is the sum of both leagues:")
    print(f"{nl_count} + {al_count} = {total_count}")

count_mlb_triple_crowns()