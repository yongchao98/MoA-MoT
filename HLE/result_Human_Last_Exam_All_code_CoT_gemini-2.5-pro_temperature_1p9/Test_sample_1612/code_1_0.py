def solve_mlb_triple_crown():
    """
    Calculates the total number of unique MLB Triple Crown winners.

    The function first defines the number of unique players who have won the
    batting Triple Crown and the number of unique players who have won the
    pitching Triple Crown. It then adds these two numbers together to find
    the total number of unique winners.
    """
    # According to historical MLB records, there are 16 unique players
    # who have won the batting Triple Crown. This includes players who
    # won it multiple times, such as Rogers Hornsby and Ted Williams,
    # who are each counted only once.
    num_batting_winners = 16

    # Similarly, there are 25 unique players who have won the pitching
    # Triple Crown. Players like Sandy Koufax and Walter Johnson, who
    # won it three times each, are only counted once.
    num_pitching_winners = 25

    # The total number of winners is the sum of unique batting winners
    # and unique pitching winners.
    total_winners = num_batting_winners + num_pitching_winners

    print("The total number of MLB Triple Crown winners is found by adding the unique batting winners and unique pitching winners.")
    print(f"Unique Batting Triple Crown Winners: {num_batting_winners}")
    print(f"Unique Pitching Triple Crown Winners: {num_pitching_winners}")
    print("\nFinal Equation:")
    print(f"{num_batting_winners} (Batters) + {num_pitching_winners} (Pitchers) = {total_winners}")

    print(f"\nThere are a total of {total_winners} unique Triple Crown winners in MLB history.")


solve_mlb_triple_crown()