# Step 1 & 2: Define the point values for each part of the meld.

# My meld points
my_run = 15 # A family (run) is worth 15 points.
my_eight_aces = 100 # Having all 8 aces is worth 100 points.

# Partner's meld points
partner_pinochle = 4 # A pinochle (Queen of Spades and Jack of Diamonds) is 4 points.
partner_nines_of_trump = 2 # Each 9 of trump is 1 point, so two are 2 points.

# Step 3: Calculate total meld points
my_total_meld = my_run + my_eight_aces
partner_total_meld = partner_pinochle + partner_nines_of_trump
total_meld = my_total_meld + partner_total_meld

# Step 4: Calculate trick points.
# With 8 aces and perfect play, my team wins all 12 tricks.
# This means we capture all "counter" cards (Aces, 10s, Kings).
# In a 48-card Pinochle deck, there are 8 of each.
# Each counter is worth 1 point.
aces_in_tricks = 8
tens_in_tricks = 8
kings_in_tricks = 8
counter_points = aces_in_tricks + tens_in_tricks + kings_in_tricks

# There is also a 1-point bonus for winning the last trick.
last_trick_bonus = 1
total_trick_points = counter_points + last_trick_bonus

# Step 5: Calculate the final total score.
final_score = total_meld + total_trick_points

# Step 6: Print the final equation as requested.
print("The total points are calculated as follows:")
print(f"{my_run} (My Run) + {my_eight_aces} (My 8 Aces) + "
      f"{partner_pinochle} (Partner's Pinochle) + {partner_nines_of_trump} (Partner's 9s of Trump) + "
      f"{aces_in_tricks} (Aces in tricks) + {tens_in_tricks} (10s in tricks) + "
      f"{kings_in_tricks} (Kings in tricks) + {last_trick_bonus} (Last Trick) = {final_score}")
