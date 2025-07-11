def count_triple_crown_winners():
    """
    This function calculates and prints the total number of unique MLB Triple Crown winners.
    It lists the batting and pitching winners separately, counts each group,
    and then sums the counts to find the total.
    """
    # List of unique players who have won the batting Triple Crown
    batting_winners = [
        "Paul Hines", "Tip O'Neill", "Nap Lajoie", "Ty Cobb", "Heinie Zimmerman",
        "Rogers Hornsby", "Chuck Klein", "Jimmie Foxx", "Lou Gehrig", "Joe Medwick",
        "Ted Williams", "Mickey Mantle", "Frank Robinson", "Carl Yastrzemski", "Miguel Cabrera"
    ]

    # List of unique players who have won the pitching Triple Crown
    pitching_winners = [
        "Tommy Bond", "Old Hoss Radbourn", "Tim Keefe", "John Clarkson", "Amos Rusie",
        "Cy Young", "Rube Waddell", "Christy Mathewson", "Walter Johnson", "Pete Alexander",
        "Dazzy Vance", "Lefty Grove", "Lefty Gomez", "Bob Feller", "Hal Newhouser",
        "Sandy Koufax", "Jake Peavy", "Clayton Kershaw", "Justin Verlander", "Shane Bieber"
    ]

    num_batting_winners = len(batting_winners)
    num_pitching_winners = len(pitching_winners)
    total_winners = num_batting_winners + num_pitching_winners

    print("List of Unique Batting Triple Crown Winners:")
    for winner in batting_winners:
        print(f"- {winner}")
    print("\n")

    print("List of Unique Pitching Triple Crown Winners:")
    for winner in pitching_winners:
        print(f"- {winner}")
    print("\n")

    print("Calculating the total number of unique Triple Crown winners:")
    print(f"Number of Batting Winners + Number of Pitching Winners = Total Winners")
    print(f"{num_batting_winners} + {num_pitching_winners} = {total_winners}")
    print("\n")
    print(f"There are a total of {total_winners} unique Triple Crown winners in MLB history.")


if __name__ == "__main__":
    count_triple_crown_winners()