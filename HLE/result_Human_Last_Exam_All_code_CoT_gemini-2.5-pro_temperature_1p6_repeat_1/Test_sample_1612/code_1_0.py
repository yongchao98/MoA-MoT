def count_mlb_triple_crown_winners():
    """
    This function calculates and prints the total number of unique MLB Triple Crown winners.
    It considers both batting and pitching Triple Crowns.
    """
    # Source: Major League Baseball historical records.
    # List of unique players who have won the batting Triple Crown.
    # Some players, like Rogers Hornsby and Ted Williams, won it multiple times but are counted once.
    batting_winners = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Heinie Zimmerman", "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx",
        "Lou Gehrig", "Joe Medwick", "Ted Williams", "Mickey Mantle",
        "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]
    num_batting_winners = len(batting_winners)

    # List of unique players who have won the pitching Triple Crown.
    # Players like Walter Johnson, Sandy Koufax, and others won multiple times but are counted once.
    pitching_winners = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Rube Waddell", "Christy Mathewson", "Walter Johnson", "Dazzy Vance",
        "Pete Alexander", "Lefty Grove", "Lefty Gomez", "Bucky Walters",
        "Bob Feller", "Hal Newhouser", "Sandy Koufax", "Dwight Gooden",
        "Roger Clemens", "Pedro Martínez", "Randy Johnson", "Johan Santana",
        "Clayton Kershaw", "Justin Verlander", "Shane Bieber"
    ]
    num_pitching_winners = len(pitching_winners)

    # Calculate the total number of unique winners.
    total_winners = num_batting_winners + num_pitching_winners

    print(f"Number of unique batting Triple Crown winners: {num_batting_winners}")
    print(f"Number of unique pitching Triple Crown winners: {num_pitching_winners}")
    print(f"Total MLB Triple Crown Winners = {num_batting_winners} + {num_pitching_winners} = {total_winners}")

count_mlb_triple_crown_winners()