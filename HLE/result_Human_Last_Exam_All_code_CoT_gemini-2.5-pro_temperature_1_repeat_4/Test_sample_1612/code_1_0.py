def count_mlb_triple_crown_winners():
    """
    Calculates and prints the total number of unique MLB Triple Crown winners.
    A Triple Crown is won by a batter leading their league in batting average, home runs, and RBIs,
    or a pitcher leading their league in wins, ERA, and strikeouts.
    """

    # List of unique players who have won the batting Triple Crown
    # Note: Rogers Hornsby and Ted Williams each won it twice.
    batting_winners = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx", "Lou Gehrig", "Joe Medwick",
        "Ted Williams", "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]

    # List of unique players who have won the pitching Triple Crown
    # Note: Several pitchers won it multiple times (e.g., Koufax, Johnson, Alexander).
    pitching_winners = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Cy Young", "Rube Waddell", "Christy Mathewson", "Walter Johnson", "Pete Alexander",
        "Dazzy Vance", "Lefty Grove", "Lefty Gomez", "Bucky Walters", "Bob Feller",
        "Hal Newhouser", "Sandy Koufax", "Steve Carlton", "Pedro Martínez", "Randy Johnson",
        "Johan Santana", "Clayton Kershaw", "Justin Verlander", "Shane Bieber"
    ]

    num_batting_winners = len(batting_winners)
    num_pitching_winners = len(pitching_winners)
    total_winners = num_batting_winners + num_pitching_winners

    print(f"There have been {num_batting_winners} unique batting Triple Crown winners.")
    print(f"There have been {num_pitching_winners} unique pitching Triple Crown winners.")
    print("\nTo find the total number of unique winners, we add these two numbers together:")
    print(f"{num_batting_winners} + {num_pitching_winners} = {total_winners}")

count_mlb_triple_crown_winners()