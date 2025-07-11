def count_mlb_triple_crowns():
    """
    Calculates and prints the total number of MLB Triple Crown winners.
    The Triple Crown can be won by either a batter or a pitcher.
    This function counts the number of times each has been achieved and sums them up.
    """

    # Data represents each time the Triple Crown was won in MLB history.
    # Source: MLB official records and Baseball-Reference.com
    batting_winners_count = 18
    pitching_winners_count = 35

    # Calculate the total
    total_winners = batting_winners_count + pitching_winners_count

    # Print the breakdown and the final equation
    print(f"Number of batting Triple Crown winners: {batting_winners_count}")
    print(f"Number of pitching Triple Crown winners: {pitching_winners_count}")
    print(f"\nTotal MLB Triple Crown Winners = {batting_winners_count} + {pitching_winners_count} = {total_winners}")

if __name__ == "__main__":
    count_mlb_triple_crowns()