# --- Point Values ---
# Meld points for the declarer (you)
declarer_run_points = 15      # A family/run (A-K-Q-J-10 of trump) is worth 15 points.
declarer_aces_points = 100    # Having all 8 aces (two of each suit) is worth 100 points.

# Meld points for the partner
partner_pinochle_points = 4   # A pinochle (Queen of Spades and Jack of Diamonds) is worth 4 points.
partner_dix_points = 2        # Each 9 of trump (dix) is worth 1 point. The partner has both, for a total of 2 points.

# Trick points
total_counter_points_in_deck = 24 # There are 24 counter cards (8 Aces, 8 Tens, 8 Kings), each worth 1 point.
last_trick_bonus = 1              # The team that wins the last trick gets a 1 point bonus.

# --- Calculations ---
# 1. Calculate the total meld points for the team.
total_declarer_meld = declarer_run_points + declarer_aces_points
total_partner_meld = partner_pinochle_points + partner_dix_points
total_meld_points = total_declarer_meld + total_partner_meld

# 2. Calculate the total trick points.
# With all 8 aces and a run in trump, perfect play guarantees winning every trick.
# This means your team will capture all counter points and the last trick bonus.
total_trick_points = total_counter_points_in_deck + last_trick_bonus

# 3. Calculate the grand total for the hand.
grand_total = total_meld_points + total_trick_points

# --- Output the Results ---
print("Here is the breakdown of your team's score for the hand:")

# Meld Point Calculation
print("\n--- Meld Points ---")
print(f"Your meld points (run + 8 aces): {declarer_run_points} + {declarer_aces_points} = {total_declarer_meld}")
print(f"Your partner's meld points (pinochle + dix): {partner_pinochle_points} + {partner_dix_points} = {total_partner_meld}")
print(f"Total Team Meld Points: {total_declarer_meld} + {total_partner_meld} = {total_meld_points}")

# Trick Point Calculation
print("\n--- Trick Points ---")
print("With all 8 aces, perfect play allows your team to win all tricks.")
print(f"This means you capture all counter card points plus the last trick bonus.")
print(f"Total Trick Points: {total_counter_points_in_deck} + {last_trick_bonus} = {total_trick_points}")

# Final Score Calculation
print("\n--- Final Score ---")
print("The total score is the sum of all meld and trick points.")
print(f"Final Equation: {total_meld_points} (Meld) + {total_trick_points} (Tricks) = {grand_total}")
print(f"\nYour team will earn a total of {grand_total} points this hand.")