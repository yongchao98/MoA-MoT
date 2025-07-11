def count_triple_crown_winners():
    """
    This function calculates and prints the total number of unique MLB Triple Crown winners.
    """
    # List of all instances of a player winning the Batting Triple Crown.
    # Some players, like Rogers Hornsby and Ted Williams, won it more than once.
    batting_winners_instances = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Heinie Zimmerman", "Rogers Hornsby", "Rogers Hornsby", "Chuck Klein",
        "Jimmie Foxx", "Lou Gehrig", "Joe Medwick", "Ted Williams", "Ted Williams",
        "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]

    # List of all instances of a player winning the Pitching Triple Crown.
    # Players like Walter Johnson, Sandy Koufax, and others won it multiple times.
    pitching_winners_instances = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Cy Young", "Rube Waddell", "Christy Mathewson", "Christy Mathewson",
        "Walter Johnson", "Pete Alexander", "Pete Alexander", "Walter Johnson",
        "Pete Alexander", "Walter Johnson", "Dazzy Vance", "Lefty Grove", "Lefty Grove",
        "Lefty Gomez", "Lefty Gomez", "Bucky Walters", "Bob Feller", "Hal Newhouser",
        "Sandy Koufax", "Sandy Koufax", "Sandy Koufax", "Dwight Gooden",
        "Roger Clemens", "Roger Clemens", "Pedro Martinez", "Randy Johnson",
        "Johan Santana", "Clayton Kershaw", "Justin Verlander", "Shane Bieber", "Gerrit Cole"
    ]

    # Use a set to get the count of unique players for each category
    unique_batting_winners = set(batting_winners_instances)
    num_unique_batting_winners = len(unique_batting_winners)

    unique_pitching_winners = set(pitching_winners_instances)
    num_unique_pitching_winners = len(unique_pitching_winners)

    # Calculate the total number of unique winners
    total_unique_winners = num_unique_batting_winners + num_unique_pitching_winners

    print("To find the total number of MLB Triple Crown winners, we count the unique winners from both batting and pitching.")
    print(f"Number of unique batting Triple Crown winners: {num_unique_batting_winners}")
    print(f"Number of unique pitching Triple Crown winners: {num_unique_pitching_winners}")
    print("\nFinal Calculation:")
    print(f"{num_unique_batting_winners} + {num_unique_pitching_winners} = {total_unique_winners}")

count_triple_crown_winners()