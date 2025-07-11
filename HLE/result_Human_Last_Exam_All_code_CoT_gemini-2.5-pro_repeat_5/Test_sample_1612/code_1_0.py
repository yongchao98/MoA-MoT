def count_mlb_triple_crown_winners():
    """
    Calculates and prints the total number of unique MLB Triple Crown winners.
    """
    # A list of all instances where a player won the batting Triple Crown.
    # Players who won multiple times are listed for each win.
    batting_winners_by_year = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Heinie Zimmerman", "Rogers Hornsby", "Rogers Hornsby", "Chuck Klein",
        "Jimmie Foxx", "Lou Gehrig", "Joe Medwick", "Ted Williams",
        "Ted Williams", "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski",
        "Miguel Cabrera"
    ]

    # A list of all instances where a player won the pitching Triple Crown.
    pitching_winners_by_year = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson",
        "Amos Rusie", "Cy Young", "Rube Waddell", "Christy Mathewson",
        "Christy Mathewson", "Walter Johnson", "Walter Johnson", "Walter Johnson",
        "Pete Alexander", "Pete Alexander", "Pete Alexander", "Hippo Vaughn",
        "Lefty Grove", "Lefty Grove", "Lefty Gomez", "Lefty Gomez",
        "Bucky Walters", "Hal Newhouser", "Sandy Koufax", "Sandy Koufax",
        "Sandy Koufax", "Steve Carlton", "Dwight Gooden", "Roger Clemens",
        "Roger Clemens", "Pedro Martinez", "Randy Johnson", "Johan Santana",
        "Clayton Kershaw", "Justin Verlander", "Shane Bieber"
    ]

    # Use a set to get the count of unique players for each category.
    # Note: Paul Hines' 1878 Triple Crown is recognized by some sources like Baseball-Reference,
    # but not officially by MLB. We will include him for a comprehensive count.
    num_unique_batting_winners = len(set(batting_winners_by_year))
    num_unique_pitching_winners = len(set(pitching_winners_by_year))

    # Calculate the total number of unique winners.
    total_unique_winners = num_unique_batting_winners + num_unique_pitching_winners

    print("Calculating the total number of unique MLB Triple Crown winners:")
    print(f"Number of unique batting Triple Crown winners: {num_unique_batting_winners}")
    print(f"Number of unique pitching Triple Crown winners: {num_unique_pitching_winners}")
    print(f"Total unique Triple Crown winners = {num_unique_batting_winners} + {num_unique_pitching_winners} = {total_unique_winners}")

count_mlb_triple_crown_winners()
<<<42>>>