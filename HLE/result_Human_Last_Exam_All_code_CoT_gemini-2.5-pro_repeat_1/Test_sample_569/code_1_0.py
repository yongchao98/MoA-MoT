# Step 1: Calculate the auction winner's (my) meld points.
my_run_meld = 15  # A, 10, K, Q, J of trump
my_aces_meld = 100  # 8 aces is "double aces"
my_total_meld = my_run_meld + my_aces_meld

# Step 2: Calculate the partner's meld points.
partner_pinochle_meld = 4  # Jack of Diamonds and Queen of Spades
partner_dix_meld = 20  # Two 9s of trump, at 10 points each
partner_total_meld = partner_pinochle_meld + partner_dix_meld

# Step 3: Calculate the total trick points.
# With the specified hand and "perfect play", the partnership wins all tricks.
# Points in deck: 4 suits * (2*11 for Aces + 2*10 for 10s + 2*4 for Kings + 2*3 for Queens + 2*2 for Jacks)
deck_card_points = 240
last_trick_bonus = 10
total_trick_points = deck_card_points + last_trick_bonus

# Step 4: Calculate the grand total score for the hand.
grand_total = my_total_meld + partner_total_meld + total_trick_points

# Print the final equation showing all components of the score.
print("The total points earned this hand are calculated as follows:")
print(f"{my_total_meld} (My Meld) + {partner_total_meld} (Partner's Meld) + {total_trick_points} (Trick Points) = {grand_total} Total Points")
