# Step 1: Calculate my meld points
my_run_points = 15  # A "family" or run in trump (A, 10, K, Q, J)
my_aces_points = 100 # "8 aces" is two of each ace, which is "100 aces"
my_total_meld = my_run_points + my_aces_points

# Step 2: Calculate my partner's meld points
partner_pinochle_points = 4  # One pinochle (J of Diamonds, Q of Spades)
# Both 9s of trump are called "dix", each is worth 10 points
partner_dix_points = 2 * 10
partner_total_meld = partner_pinochle_points + partner_dix_points

# Step 3: Calculate the points from winning all the tricks
# Points per suit: 2*A(11) + 2*10(10) + 2*K(4) + 2*Q(3) + 2*J(2) = 22+20+8+6+4 = 60
# Total points in 4 suits = 60 * 4 = 240
deck_card_points = 240
# Bonus for winning the last trick
last_trick_bonus = 10
# With perfect play and the overwhelming trump strength described,
# it is assumed all tricks are won.
total_trick_points = deck_card_points + last_trick_bonus

# Step 4: Calculate the grand total score for the hand
total_score = my_total_meld + partner_total_meld + total_trick_points

# Final output showing the breakdown of all points
print("The total points earned this hand are calculated as follows:")
print(f"My Meld: Run ({my_run_points}) + 100 Aces ({my_aces_points}) = {my_total_meld} points")
print(f"Partner's Meld: Pinochle ({partner_pinochle_points}) + Two Dix ({partner_dix_points}) = {partner_total_meld} points")
print(f"Trick Points: Card Points ({deck_card_points}) + Last Trick Bonus ({last_trick_bonus}) = {total_trick_points} points")
print("-" * 30)
# The final equation as requested
print(f"Final Equation: {my_run_points} (Run) + {my_aces_points} (100 Aces) + {partner_pinochle_points} (Pinochle) + {partner_dix_points} (Two Dix) + {deck_card_points} (Cards) + {last_trick_bonus} (Last Trick) = {total_score}")
print(f"Total Score = {total_score} points")
