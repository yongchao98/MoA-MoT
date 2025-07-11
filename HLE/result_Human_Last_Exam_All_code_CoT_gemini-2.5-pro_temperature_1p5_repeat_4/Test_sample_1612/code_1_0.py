# A function to solve the problem
def count_triple_crown_winners():
    """
    Calculates and prints the total number of MLB Triple Crown winners.
    This is done by counting the unique winners for both batting and pitching
    and summing them up.
    """
    # Using a set to store names ensures that players who won multiple times are only counted once.
    # Data is based on MLB official records.
    batting_winners = {
        "Tip O'Neill", "Hugh Duffy", "Nap Lajoie", "Ty Cobb", "Heinie Zimmerman",
        "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx", "Lou Gehrig", "Joe Medwick",
        "Ted Williams", "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    }

    pitching_winners = {
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Cy Young", "Rube Waddell", "Christy Mathewson", "Walter Johnson", "Pete Alexander",
        "Hippo Vaughn", "Dazzy Vance", "Lefty Grove", "Lefty Gomez", "Bucky Walters",
        "Bob Feller", "Hal Newhouser", "Sandy Koufax", "Dwight Gooden", "Roger Clemens",
        "Pedro Martinez", "Randy Johnson", "Johan Santana", "Clayton Kershaw", "Justin Verlander",
        "Shane Bieber"
    }

    num_batting_winners = len(batting_winners)
    num_pitching_winners = len(pitching_winners)
    total_winners = num_batting_winners + num_pitching_winners

    print(f"Number of unique batting Triple Crown winners: {num_batting_winners}")
    print(f"Number of unique pitching Triple Crown winners: {num_pitching_winners}")
    print("\nTotal number of Triple Crown winners in MLB history:")
    # The final print statement shows the equation with each number
    print(f"{num_batting_winners} (batting) + {num_pitching_winners} (pitching) = {total_winners}")

# Execute the function
count_triple_crown_winners()
<<<41>>>