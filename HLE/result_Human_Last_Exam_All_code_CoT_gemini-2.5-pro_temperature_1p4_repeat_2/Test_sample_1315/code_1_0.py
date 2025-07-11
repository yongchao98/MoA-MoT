def solve_tichu_score_difference():
    """
    Calculates the maximal possible score difference (X-Y) in a single round of Tichu
    where the winning team does not finish 1st and 2nd.
    """

    # --- 1. Card Points Calculation ---
    # The total value of all point cards (Kings, 10s, 5s, Dragon, Phoenix) is 100.
    # To maximize the winning team's (Team A) score advantage from cards, a player
    # from the losing team (Team B) must finish last. This transfers their card points to Team A.
    # The optimal finishing order is: 1st (Team A), 2nd (Team B), 3rd (Team A), 4th (Team B).
    #
    # To maximize the score difference, Team A must capture all positive point cards
    # (worth 100 points) and ensure Team B captures only the negative point card (the Phoenix, worth -25).
    # This happens if the Team B player who finished 2nd is the one to win the trick with the Phoenix.
    winning_team_card_points = 125
    losing_team_card_points = -25

    # --- 2. Tichu Call Points Calculation ---
    # To maximize the difference from calls, the winning team must make a successful high-value call,
    # and the losing team must make failed calls.
    #
    # The player from Team A who went out first makes a successful "Grand Tichu" call.
    winning_team_tichu_points = 200
    # Both players on Team B make "Tichu" calls. Since they do not go out first, they both fail.
    losing_team_tichu_points = -100 + -100

    # --- 3. Total Score and Difference Calculation ---
    # The winning team's score, X, is the sum of their card and tichu points.
    X = winning_team_card_points + winning_team_tichu_points
    # The losing team's score, Y, is the sum of their card and tichu points.
    Y = losing_team_card_points + losing_team_tichu_points

    # The final difference, X - Y.
    X_minus_Y = X - Y

    # --- Printing the step-by-step calculation ---
    print("This script calculates the maximal score difference (X-Y) in a game of Tichu.")
    print("X is the winning team's score, and Y is the losing team's score.")
    print("The winning team is defined as the team of the player who went out first.\n")

    print(f"Step 1: Calculate the maximum card point difference.")
    print(f" - The winning team scores: {winning_team_card_points} points.")
    print(f" - The losing team scores: {losing_team_card_points} points.")
    print(f" - Card Point Difference = {winning_team_card_points} - ({losing_team_card_points}) = {winning_team_card_points - losing_team_card_points}\n")
    
    print(f"Step 2: Calculate the maximum Tichu call point difference.")
    print(f" - The winning team makes a successful Grand Tichu: +{winning_team_tichu_points} points.")
    print(f" - The losing team has two failed Tichu calls: -100 + -100 = {losing_team_tichu_points} points.")
    print(f" - Tichu Point Difference = {winning_team_tichu_points} - ({losing_team_tichu_points}) = {winning_team_tichu_points - losing_team_tichu_points}\n")

    print("Step 3: Calculate the final scores X and Y.")
    print(f" - Winning Team Score (X) = Card Points + Tichu Points")
    print(f"   X = {winning_team_card_points} + {winning_team_tichu_points} = {X}")
    print(f" - Losing Team Score (Y) = Card Points + Tichu Points")
    print(f"   Y = {losing_team_card_points} + ({losing_team_tichu_points}) = {Y}\n")

    print("Final Result: The maximal possible value of X - Y is:")
    print(f"   X - Y = {X} - ({Y}) = {X_minus_Y}")

solve_tichu_score_difference()
<<<550>>>