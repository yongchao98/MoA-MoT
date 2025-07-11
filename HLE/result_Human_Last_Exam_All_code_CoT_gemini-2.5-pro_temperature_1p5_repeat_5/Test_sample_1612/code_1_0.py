def count_triple_crown_winners():
    """
    Calculates and prints the total number of unique MLB Triple Crown winners.
    The Triple Crown can be won by either a batter or a pitcher. This function
    sums the unique winners from both categories.
    """
    # Number of unique players who have won the batting Triple Crown
    # (leading their league in batting average, home runs, and RBIs).
    # Some players, like Rogers Hornsby, won it multiple times but are counted once.
    unique_batting_winners = 15

    # Number of unique players who have won the pitching Triple Crown
    # (leading their league in wins, strikeouts, and ERA).
    # Players like Walter Johnson and Sandy Koufax won multiple times but are counted once.
    unique_pitching_winners = 23

    # Total number of unique players who have won a Triple Crown.
    # Since no player has won both a batting and a pitching Triple Crown,
    # we can simply add the two numbers together.
    total_unique_winners = unique_batting_winners + unique_pitching_winners

    print("To find the total number of MLB Triple Crown winners, we add the unique batting winners and unique pitching winners.")
    print(f"Unique Batting Triple Crown Winners: {unique_batting_winners}")
    print(f"Unique Pitching Triple Crown Winners: {unique_pitching_winners}")
    print("The total number of unique Triple Crown winners is the sum of these two groups:")
    print(f"{unique_batting_winners} + {unique_pitching_winners} = {total_unique_winners}")

count_triple_crown_winners()