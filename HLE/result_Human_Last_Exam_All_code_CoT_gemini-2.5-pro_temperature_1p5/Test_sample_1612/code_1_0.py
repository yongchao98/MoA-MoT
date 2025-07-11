def count_triple_crown_winners():
    """
    Calculates and prints the total number of unique MLB Triple Crown winners.
    """
    # Batting Triple Crown winners. Some players, like Rogers Hornsby and Ted Williams, won more than once.
    batting_winners = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Heinie Zimmerman", "Rogers Hornsby", "Rogers Hornsby", "Jimmie Foxx",
        "Chuck Klein", "Lou Gehrig", "Joe Medwick", "Ted Williams", "Ted Williams",
        "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]

    # Pitching Triple Crown winners. Many pitchers have won multiple times.
    pitching_winners = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Cy Young", "Christy Mathewson", "Rube Waddell", "Christy Mathewson",
        "Walter Johnson", "Grover Cleveland Alexander", "Grover Cleveland Alexander",
        "Walter Johnson", "Lefty Grove", "Lefty Grove", "Grover Cleveland Alexander",
        "Dazzy Vance", "Walter Johnson", "Lefty Gomez", "Lefty Gomez", "Bucky Walters",
        "Hal Newhouser", "Sandy Koufax", "Sandy Koufax", "Sandy Koufax",
        "Jake Peavy", "Clayton Kershaw", "Justin Verlander", "Shane Bieber"
    ]

    # Use a set to get the count of unique players for each category
    unique_batting_winners = set(batting_winners)
    unique_pitching_winners = set(pitching_winners)

    num_batting_winners = len(unique_batting_winners)
    num_pitching_winners = len(unique_pitching_winners)

    # The total number of unique winners is the sum of unique winners from both categories.
    # No player has ever won both a batting and a pitching Triple Crown.
    total_unique_winners = num_batting_winners + num_pitching_winners

    print(f"Number of unique batting Triple Crown winners: {num_batting_winners}")
    print(f"Number of unique pitching Triple Crown winners: {num_pitching_winners}")
    print("\nTotal number of unique MLB Triple Crown winners:")
    # The final equation is printed as requested
    print(f"{num_batting_winners} + {num_pitching_winners} = {total_unique_winners}")

count_triple_crown_winners()