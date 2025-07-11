def solve_pinochle_hand():
    """
    Calculates the total score for a given Pinochle hand scenario.
    """

    # 1. Meld Points Calculation
    # Your meld (auction winner)
    my_run_meld = 15  # A run (A, 10, K, Q, J of trump)
    my_aces_meld = 100 # 8 aces (two of each suit)
    my_total_meld = my_run_meld + my_aces_meld

    # Partner's meld
    partner_pinochle_meld = 4 # One pinochle (Q of Spades, J of Diamonds)
    partner_nines_meld = 2 * 10 # Both 9s of trump (each is a 'dix' worth 10)
    partner_total_meld = partner_pinochle_meld + partner_nines_meld

    # Total team meld points
    total_meld_points = my_total_meld + partner_total_meld

    # 2. Trick Points Calculation
    # Point values for cards that score in tricks
    ace_value = 11
    ten_value = 10
    king_value = 4
    queen_value = 3
    jack_value = 2

    # Total points in a single suit set (2 of each card)
    # (11+10+4+3+2) * 2 = 60 points per suit
    points_per_suit = 2 * (ace_value + ten_value + king_value + queen_value + jack_value)

    # Total points in the 48-card deck (4 suits)
    # 60 * 4 = 240
    total_trick_points = points_per_suit * 4

    # 3. Bonus Points Calculation
    # Bonus for winning the last trick
    last_trick_bonus = 10

    # 4. Final Score Calculation
    # Under the assumption of perfect play with a dominant hand,
    # your team wins all tricks and thus all trick points plus the bonus.
    total_score = total_meld_points + total_trick_points + last_trick_bonus

    # Output the breakdown of the calculation
    print("Calculating the total score for the hand:")
    print("-" * 40)
    print("1. Total Meld Points")
    print(f"   - Your Meld (Run + 8 Aces): {my_run_meld} + {my_aces_meld} = {my_total_meld}")
    print(f"   - Partner's Meld (Pinochle + 2 Nines): {partner_pinochle_meld} + {partner_nines_meld} = {partner_total_meld}")
    print(f"   - Subtotal Meld: {my_total_meld} + {partner_total_meld} = {total_meld_points} points\n")

    print("2. Total Points from Play")
    print(f"   - Points in all cards in the deck: {total_trick_points}")
    print(f"   - Bonus for winning the last trick: {last_trick_bonus}")
    print(f"   - Subtotal from Play: {total_trick_points} + {last_trick_bonus} = {total_trick_points + last_trick_bonus} points\n")
    
    print("Final Score Calculation:")
    print(f"   {total_meld_points} (Meld) + {total_trick_points} (Tricks) + {last_trick_bonus} (Bonus) = {total_score} points")

    # Final answer in the required format
    print(f"\n<<<{total_score}>>>")

solve_pinochle_hand()