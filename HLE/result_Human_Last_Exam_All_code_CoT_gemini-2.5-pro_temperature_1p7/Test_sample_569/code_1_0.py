import sys

def solve_pinochle_score():
    """
    This function calculates the total score for a Pinochle hand based on the user's description.
    """
    # Step 1: Calculate meld points for the team.
    
    # My meld points
    # A "family" or run in trump is 15 points.
    my_run_points = 15
    # Having all 8 aces ("100 aces") is 100 points.
    my_aces_points = 100
    my_total_meld = my_run_points + my_aces_points
    
    # Partner's meld points
    # A "pinochle" (Q of spades, J of diamonds) is 4 points.
    partner_pinochle_points = 4
    # Each 9 of trump ("dix") is worth 10 points. Both 9s are 20 points.
    partner_dix_points = 2 * 10
    partner_total_meld = partner_pinochle_points + partner_dix_points
    
    # Total team meld points
    total_meld_points = my_total_meld + partner_total_meld
    
    # Step 2: Calculate trick-taking points.
    
    # With your run, all 8 aces, and your partner's two 9s of trump,
    # your team holds all 8 trump cards. With perfect play, you win all 12 tricks.
    # Therefore, you collect all card points and the last trick bonus.
    
    # Point values for cards
    ace_value = 11
    ten_value = 10
    king_value = 4
    queen_value = 3
    jack_value = 2
    nine_value = 0
    
    # There are 8 of each rank in a 48-card Pinochle deck.
    cards_per_rank = 8
    
    total_card_points = cards_per_rank * (ace_value + ten_value + king_value + queen_value + jack_value + nine_value)
    
    # Bonus for winning the last trick
    last_trick_bonus = 10
    
    # Total points from trick-taking phase
    total_trick_points = total_card_points + last_trick_bonus
    
    # Step 3: Calculate the final total score.
    final_score = total_meld_points + total_trick_points
    
    # Output the final calculation as an equation, showing each number.
    print("Your team's score is the sum of your meld, your partner's meld, and the points from winning all the tricks.")
    print(f"Final Equation:")
    print(f"{my_total_meld} (Your Meld) + {partner_total_meld} (Partner's Meld) + {total_card_points} (Card Points) + {last_trick_bonus} (Last Trick Bonus) = {final_score}")
    
    # Also output the final answer in the required format to a file.
    # The actual output for the user is handled by the print statements above.
    sys.stdout = open('output.txt', 'w')
    print(f'<<<{final_score}>>>')
    sys.stdout.close()
    sys.stdout = sys.__stdout__

solve_pinochle_score()