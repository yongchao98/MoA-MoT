# Point values for various Pinochle melds and tricks

# My meld points
my_run_points = 15
my_aces_points = 100

# Partner's meld points
partner_pinochle_points = 4
# A "dix" (9 of trump) is 10 points. The partner has two.
partner_dix_points = 10 * 2

# Trick points
# In a 48-card deck, there are 24 point-cards (A, 10, K).
# Winning the last trick gives an extra point.
# Assuming perfect play with this hand means winning all tricks.
total_trick_points = 24 + 1

# Calculate the total score
total_points = my_run_points + my_aces_points + partner_pinochle_points + partner_dix_points + total_trick_points

# Print the final equation and the result
print("Calculating the total score for the hand:")
print(f"My Meld (Run + 1000 Aces): {my_run_points} + {my_aces_points}")
print(f"Partner's Meld (Pinochle + 2 Dix): {partner_pinochle_points} + {partner_dix_points}")
print(f"Total Trick Points (All tricks taken): {total_trick_points}")
print("---")
print(f"Total Score = {my_run_points} + {my_aces_points} + {partner_pinochle_points} + {partner_dix_points} + {total_trick_points} = {total_points}")
