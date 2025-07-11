# Step 1: Calculate my meld points
my_run_points = 15
my_double_aces_points = 100
my_total_meld = my_run_points + my_double_aces_points

# Step 2: Calculate my partner's meld points
partner_pinochle_points = 4
partner_dix_points = 2 # 1 point for each of the two 9s of trump
partner_total_meld = partner_pinochle_points + partner_dix_points

# Step 3: Calculate the trick points
# In a 48-card deck, there are 8 Aces, 8 Tens, and 8 Kings.
# Each is worth 1 point, for 24 total card points.
# We assume we win all tricks due to "perfect play".
trick_card_points = 24
last_trick_bonus = 2
total_trick_points = trick_card_points + last_trick_bonus

# Step 4: Calculate the final total score
final_score = my_total_meld + partner_total_meld + total_trick_points

# Print the breakdown and the final equation
print(f"My meld points: {my_run_points} (run) + {my_double_aces_points} (double aces) = {my_total_meld}")
print(f"Partner's meld points: {partner_pinochle_points} (pinochle) + {partner_dix_points} (dix) = {partner_total_meld}")
print(f"Trick points: {trick_card_points} (counters) + {last_trick_bonus} (last trick) = {total_trick_points}")
print("-" * 30)
print("Final Score Equation:")
print(f"{my_total_meld} + {partner_total_meld} + {total_trick_points} = {final_score}")
