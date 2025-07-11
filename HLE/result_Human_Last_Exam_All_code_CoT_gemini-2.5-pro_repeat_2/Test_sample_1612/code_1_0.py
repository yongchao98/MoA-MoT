def solve_mlb_triple_crown():
    """
    Calculates and prints the total number of unique MLB Triple Crown winners.

    The Triple Crown can be won by a batter (leading the league in batting average,
    home runs, and RBIs) or a pitcher (leading in wins, strikeouts, and ERA).
    This function sums the unique winners from both categories.
    """

    # There are 16 unique players who have won the batting Triple Crown.
    num_batting_winners = 16

    # There are 22 unique players who have won the pitching Triple Crown.
    num_pitching_winners = 22

    # Calculate the total number of unique winners.
    total_winners = num_batting_winners + num_pitching_winners

    print("To find the total number of unique Triple Crown winners in MLB, we add the number of batting winners and pitching winners.")
    print("\nNumber of unique batting Triple Crown winners: {}".format(num_batting_winners))
    print("Number of unique pitching Triple Crown winners: {}".format(num_pitching_winners))
    print("\nFinal calculation:")
    print("{} + {} = {}".format(num_batting_winners, num_pitching_winners, total_winners))

solve_mlb_triple_crown()
<<<38>>>