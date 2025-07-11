# Step 1: Define the values for meld points.
my_run = 15
my_aces = 100
partner_pinochle = 4
partner_nines = 10

# Calculate total meld points
total_meld = my_run + my_aces + partner_pinochle + partner_nines

# Step 2: Define points from tricks.
# In a 48-card deck, there are 8 Aces, 8 Tens, and 8 Kings.
# Each is worth 1 point, for a total of 24 trick points.
# Having all 8 aces and playing perfectly guarantees winning all these points.
trick_points = 24

# Step 3: Define the bonus for winning the last trick.
# This is guaranteed with the hand described.
last_trick_bonus = 1

# Step 4: Calculate the total score for the hand.
total_score = total_meld + trick_points + last_trick_bonus

# Print the breakdown of the score calculation.
print("Calculating the total score:")
print(f"Your Meld (Run + 8 Aces): {my_run} + {my_aces} = {my_run + my_aces}")
print(f"Partner's Meld (Pinochle + 2 Nines): {partner_pinochle} + {partner_nines} = {partner_pinochle + partner_nines}")
print(f"Total Meld Points: {total_meld}")
print(f"Points from Tricks: {trick_points}")
print(f"Bonus for Last Trick: {last_trick_bonus}")
print("-" * 20)
print(f"Final Equation: {my_run} (run) + {my_aces} (aces) + {partner_pinochle} (pinochle) + {partner_nines} (nines) + {trick_points} (tricks) + {last_trick_bonus} (last trick)")
print(f"Total Points Earned: {total_score}")

<<<154>>>