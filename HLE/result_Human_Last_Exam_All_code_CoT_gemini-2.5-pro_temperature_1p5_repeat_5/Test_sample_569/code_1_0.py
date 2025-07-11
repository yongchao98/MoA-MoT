# Calculate your meld points
my_run_points = 15
my_aces_around_points = 100
my_total_meld = my_run_points + my_aces_around_points

# Calculate your partner's meld points
partner_pinochle_points = 4
partner_dix_points = 10 * 2
partner_total_meld = partner_pinochle_points + partner_dix_points

# Calculate trick points
# There are 8 Aces, 8 Tens, and 8 Kings, each worth 1 point in tricks.
# 1 additional point is awarded for the last trick.
# Playing a perfect hand with this meld means taking all tricks.
trick_points = (8 + 8 + 8) + 1

# Calculate the total score for the hand
total_score = my_total_meld + partner_total_meld + trick_points

# Print the breakdown of the score
print("Points Breakdown:")
print(f"Your Meld (Run + Aces): {my_run_points} + {my_aces_around_points} = {my_total_meld}")
print(f"Partner's Meld (Pinochle + Dixes): {partner_pinochle_points} + {partner_dix_points} = {partner_total_meld}")
print(f"Trick Points (Counters + Last Trick): {trick_points}")
print("-" * 20)
print("Final Equation:")
print(f"({my_run_points} + {my_aces_around_points}) + ({partner_pinochle_points} + {partner_dix_points}) + {trick_points} = {total_score}")
print(f"Total Score for the hand: {total_score}")
