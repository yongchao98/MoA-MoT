# Pinochle score calculation

# Step 1 & 2: Define the meld point values for each hand
my_run_meld = 15
my_aces_meld = 100
partner_pinochle_meld = 4
partner_dix_meld = 2 # 1 point for each of the two 9s of trump

# Step 3: Calculate total meld points
my_total_meld = my_run_meld + my_aces_meld
partner_total_meld = partner_pinochle_meld + partner_dix_meld
total_meld_points = my_total_meld + partner_total_meld

# Step 4: Calculate total trick points
# With perfect play, the team wins all tricks.
# There are 24 counter cards (8 Aces, 8 Tens, 8 Kings), each worth 1 point.
trick_card_points = 24
# There is a 1 point bonus for winning the last trick.
last_trick_bonus = 1
total_trick_points = trick_card_points + last_trick_bonus

# Step 5: Calculate the final total score
final_score = total_meld_points + total_trick_points

# Print the final breakdown as a single equation
print("The total points earned are calculated by summing the components of the team's score.")
print("The final equation is:")
print(f"{my_run_meld} (My Run) + {my_aces_meld} (My 8 Aces) + {partner_pinochle_meld} (Partner's Pinochle) + {partner_dix_meld} (Partner's 2 Nines) + {trick_card_points} (Trick Counters) + {last_trick_bonus} (Last Trick Bonus) = {final_score}")
print(f"\nThe total score for the hand is: {final_score}")
