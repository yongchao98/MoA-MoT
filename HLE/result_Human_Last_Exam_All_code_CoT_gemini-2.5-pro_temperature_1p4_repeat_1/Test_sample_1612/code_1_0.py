def count_triple_crown_winners():
    """
    This function calculates and displays the total number of unique MLB
    players who have won either a batting or a pitching Triple Crown.
    """

    # List of unique players who have won the Batting Triple Crown.
    # Note: Players who won multiple times (e.g., Rogers Hornsby, Ted Williams) are listed once.
    batting_winners = [
        "Paul Hines", "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb",
        "Heinie Zimmerman", "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx",
        "Lou Gehrig", "Joe Medwick", "Ted Williams", "Mickey Mantle",
        "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]

    # List of unique players who have won the Pitching Triple Crown (modern era, since 1900).
    # Note: Players who won multiple times (e.g., Walter Johnson, Sandy Koufax) are listed once.
    pitching_winners = [
        "Christy Mathewson", "Cy Young", "Rube Waddell", "Walter Johnson",
        "Pete Alexander", "Lefty Grove", "Lefty Gomez", "Bucky Walters",
        "Hal Newhouser", "Sandy Koufax", "Dwight Gooden", "Roger Clemens",
        "Pedro Martinez", "Randy Johnson", "Johan Santana", "Clayton Kershaw",
        "Justin Verlander", "Shane Bieber", "Gerrit Cole"
    ]

    num_batting_winners = len(batting_winners)
    num_pitching_winners = len(pitching_winners)

    # In MLB history, no player has won both a batting and a pitching Triple Crown,
    # so the total number of unique winners is the sum of the two lists.
    total_unique_winners = num_batting_winners + num_pitching_winners

    print("--- MLB Triple Crown Winners ---")
    print(f"Number of unique Batting Triple Crown Winners: {num_batting_winners}")
    print(f"Number of unique Pitching Triple Crown Winners: {num_pitching_winners}")
    print("\n--- Final Calculation ---")
    print("The total number of unique players who have won a Triple Crown is:")
    # The final equation showing each number as requested
    print(f"{num_batting_winners} (Batting) + {num_pitching_winners} (Pitching) = {total_unique_winners}")

count_triple_crown_winners()
<<<35>>>