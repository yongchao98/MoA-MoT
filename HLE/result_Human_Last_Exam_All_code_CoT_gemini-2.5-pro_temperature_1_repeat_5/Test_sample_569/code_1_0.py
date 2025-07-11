def solve_pinochle_hand():
    """
    Calculates and prints the total score for a specific Pinochle hand.
    """
    # Step 1 & 2: Calculate meld points for both players.
    my_run_meld = 15
    my_aces_meld = 100
    my_total_meld = my_run_meld + my_aces_meld

    partner_pinochle_meld = 4
    partner_nines_meld = 20  # 10 points for each of the two 9s of trump
    partner_total_meld = partner_pinochle_meld + partner_nines_meld

    # Step 3: Calculate total team meld points.
    team_meld_total = my_total_meld + partner_total_meld

    # Step 4: Calculate trick points.
    # With 8 aces, perfect play ensures winning all tricks.
    # Total points in a 48-card deck:
    # 8 Aces (11 pts), 8 Tens (10 pts), 8 Kings (4 pts), 8 Queens (3 pts), 8 Jacks (2 pts)
    # (8*11) + (8*10) + (8*4) + (8*3) + (8*2) = 88 + 80 + 32 + 24 + 16 = 240
    total_card_points_in_deck = 240
    last_trick_bonus = 10
    team_trick_points = total_card_points_in_deck + last_trick_bonus

    # Step 5: Calculate the final score.
    final_score = team_meld_total + team_trick_points

    # Print the detailed breakdown of the score calculation.
    print("Calculating the total points for the hand:")
    print(f"Our team's meld points:")
    print(f"  My meld: {my_run_meld} (run) + {my_aces_meld} (8 aces) = {my_total_meld} points")
    print(f"  Partner's meld: {partner_pinochle_meld} (pinochle) + {partner_nines_meld} (two 9s of trump) = {partner_total_meld} points")
    print(f"  Total Meld Points: {my_total_meld} + {partner_total_meld} = {team_meld_total} points")
    print("\nOur team's trick points:")
    print("  With perfect play, our team wins all tricks, securing all card points and the last trick bonus.")
    print(f"  Total Trick Points: {total_card_points_in_deck} (all card points) + {last_trick_bonus} (last trick bonus) = {team_trick_points} points")
    print("\nFinal Score Calculation:")
    print(f"Total Points = Total Meld Points + Total Trick Points")
    print(f"Total Points = {team_meld_total} + {team_trick_points} = {final_score}")

solve_pinochle_hand()
<<<389>>>