# Plan: Calculate the total points by summing meld points and trick points.

# 1. Calculate Meld Points
# My meld points
my_run_points = 15
my_eight_aces_points = 100

# Partner's meld points
partner_pinochle_points = 4
partner_two_nines_points = 20 # 10 points for each 9 of trump

# 2. Calculate Trick Points
# The total points available from countable cards in a 48-card deck is 240.
# The bonus for winning the last trick is 10.
# With a run, 8 aces, and perfect play, our team wins all tricks.
total_trick_points = 240
last_trick_bonus = 10

# 3. Calculate Final Score
# Sum all the points together.
total_score = (my_run_points + 
               my_eight_aces_points + 
               partner_pinochle_points + 
               partner_two_nines_points + 
               total_trick_points + 
               last_trick_bonus)

# Print the final equation with each component value
print("The total points earned are calculated by summing the meld and trick points.")
print("Final score equation:")
print(f"{my_run_points} (my run) + {my_eight_aces_points} (my 8 aces) + {partner_pinochle_points} (partner's pinochle) + {partner_two_nines_points} (partner's two nines) + {total_trick_points} (card points) + {last_trick_bonus} (last trick) = {total_score}")
