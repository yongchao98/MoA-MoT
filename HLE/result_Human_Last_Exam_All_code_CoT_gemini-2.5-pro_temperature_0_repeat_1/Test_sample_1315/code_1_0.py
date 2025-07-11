def solve_tichu_puzzle():
    """
    Calculates the maximal possible value of X-Y in a single Tichu round
    under the given conditions.
    """

    # Step 1: Maximize the score difference from Tichu calls.
    # The winning team (A) makes a successful Grand Tichu.
    winning_team_tichu_bonus = 200
    # The losing team (B) makes a failed Grand Tichu.
    losing_team_tichu_bonus = -200
    # The maximum possible difference from Tichu calls.
    tichu_call_difference = winning_team_tichu_bonus - losing_team_tichu_bonus

    # Step 2: Maximize the score difference from card play.
    # This is maximized when a player from the losing team finishes last,
    # and their partner (the other player on the losing team) gets the
    # minimum possible score from tricks.
    # The minimum score is from capturing the Phoenix (-25) and no other points.
    min_player_card_score = -25
    # The formula for the card play difference in this optimal scenario is:
    # 100 - 2 * (score of the non-last player on the losing team)
    card_play_difference = 100 - 2 * min_player_card_score

    # Step 3: Calculate the total maximal difference.
    # This is the sum of the difference from card play and Tichu calls.
    total_max_difference = card_play_difference + tichu_call_difference

    # Print the step-by-step calculation.
    print("The maximal possible value of X-Y is calculated by combining the maximum possible score differences from both Tichu calls and card play.")
    print("-" * 30)

    print("1. Difference from Tichu Calls:")
    print(f"   - Winning team's successful Grand Tichu bonus: {winning_team_tichu_bonus}")
    print(f"   - Losing team's failed Grand Tichu penalty: {losing_team_tichu_bonus}")
    print(f"   - Maximum difference from calls = {winning_team_tichu_bonus} - ({losing_team_tichu_bonus}) = {tichu_call_difference}")
    print("-" * 30)

    print("2. Difference from Card Play:")
    print("   This is maximized when a player from the losing team finishes last.")
    print(f"   Their partner must get the minimum possible score, which is {min_player_card_score} points.")
    print(f"   - Maximum difference from card play = 100 - 2 * ({min_player_card_score}) = {card_play_difference}")
    print("-" * 30)

    print("3. Total Maximal Difference (X - Y):")
    print(f"   Total Difference = (Card Play Difference) + (Tichu Call Difference)")
    print(f"   Total Difference = {card_play_difference} + {tichu_call_difference} = {total_max_difference}")

solve_tichu_puzzle()