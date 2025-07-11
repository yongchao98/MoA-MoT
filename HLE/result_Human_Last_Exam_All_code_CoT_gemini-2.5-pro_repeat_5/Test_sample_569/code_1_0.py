# --- Meld Point Calculation ---

# My meld points
my_run_points = 15
my_aces_points = 100

# Partner's meld points
partner_pinochle_points = 4
partner_nines_points = 2

# --- Trick Point Calculation ---

# With perfect play, your team wins all tricks.
# Total points in the deck (240) + bonus for the last trick (10).
total_trick_points = 250

# --- Final Score Calculation ---

# Sum of all points
total_score = my_run_points + my_aces_points + partner_pinochle_points + partner_nines_points + total_trick_points

# Print the final equation showing how the score is calculated.
# Each component of the score is explicitly shown in the sum.
print("The total score is calculated by summing the meld and trick points:")
print(f"{my_run_points} (my run) + {my_aces_points} (my 8 aces) + {partner_pinochle_points} (partner's pinochle) + {partner_nines_points} (partner's nines) + {total_trick_points} (all trick points) = {total_score}")