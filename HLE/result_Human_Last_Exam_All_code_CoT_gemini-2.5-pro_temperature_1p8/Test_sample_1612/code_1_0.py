def count_triple_crown_winners():
    """
    This function calculates and prints the total number of unique MLB Triple Crown winners.
    It considers both the batting and pitching Triple Crowns.
    """
    # A player wins the batting Triple Crown by leading the league in
    # batting average, home runs, and runs batted in (RBI).
    # This list contains all instances of a batting Triple Crown being won.
    batting_winners_by_instance = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Rogers Hornsby", "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx",
        "Lou Gehrig", "Joe Medwick", "Ted Williams", "Ted Williams",
        "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]

    # A player wins the pitching Triple Crown by leading the league in
    # wins, strikeouts, and earned run average (ERA).
    # This list contains all instances of a pitching Triple Crown being won.
    pitching_winners_by_instance = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Cy Young", "Rube Waddell", "Christy Mathewson", "Christy Mathewson",
        "Walter Johnson", "Walter Johnson", "Walter Johnson", "Pete Alexander",
        "Pete Alexander", "Pete Alexander", "Hippo Vaughn", "Dazzy Vance",
        "Lefty Grove", "Lefty Grove", "Bucky Walters", "Lefty Gomez", "Hal Newhouser",
        "Sandy Koufax", "Sandy Koufax", "Sandy Koufax", "Dwight Gooden",
        "Roger Clemens", "Roger Clemens", "Pedro Martínez", "Randy Johnson",
        "Johan Santana", "Clayton Kershaw", "Justin Verlander", "Shane Bieber",
        "Gerrit Cole"
    ]

    # Use a set to get the count of unique players for each category
    unique_batting_winners = set(batting_winners_by_instance)
    unique_pitching_winners = set(pitching_winners_by_instance)

    num_unique_batting = len(unique_batting_winners)
    num_unique_pitching = len(unique_pitching_winners)

    # Combine the sets to find the total number of unique players.
    # The union operator '|' automatically handles duplicates if any player won both.
    # (No player in history has won both a batting and pitching Triple Crown).
    total_unique_winners = len(unique_batting_winners | unique_pitching_winners)

    print(f"There are {num_unique_batting} unique batting Triple Crown winners and {num_unique_pitching} unique pitching Triple Crown winners.")
    print(f"The total number of unique Triple Crown winners in MLB history is {num_unique_batting} + {num_unique_pitching} = {total_unique_winners}.")

# Run the function to display the answer.
count_triple_crown_winners()
<<<41>>>