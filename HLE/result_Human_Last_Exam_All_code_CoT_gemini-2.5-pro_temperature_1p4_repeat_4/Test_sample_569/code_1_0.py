# 1. Calculate meld points for each player and the team.
my_run_meld = 15
my_aces_meld = 100 # For 8 aces (double aces around)
my_total_meld = my_run_meld + my_aces_meld

partner_pinochle_meld = 4
partner_nines_meld = 2 # For both 9s of trump
partner_total_meld = partner_pinochle_meld + partner_nines_meld

total_meld_points = my_total_meld + partner_total_meld

# 2. Calculate points from tricks.
# With all 8 aces and top trump, perfect play means winning all tricks.
# There are 8 Aces, 8 Tens, and 8 Kings, each worth 1 point.
total_trick_points = 24

# 3. Bonus for the last trick.
last_trick_bonus = 1

# 4. Calculate the total score for the hand.
total_score = total_meld_points + total_trick_points + last_trick_bonus

print("Calculating the total score for the hand:")
print("-" * 40)
print(f"My Meld: {my_run_meld} (Run) + {my_aces_meld} (8 Aces) = {my_total_meld} points")
print(f"Partner's Meld: {partner_pinochle_meld} (Pinochle) + {partner_nines_meld} (Two 9s) = {partner_total_meld} points")
print(f"Total Team Meld Points: {total_meld_points}")
print("-" * 40)
print(f"Points from winning all tricks (counters): {total_trick_points}")
print(f"Bonus for winning the last trick: {last_trick_bonus}")
print("-" * 40)
print("Final Score Equation:")
print(f"{total_meld_points} (meld) + {total_trick_points} (tricks) + {last_trick_bonus} (last trick bonus) = {total_score}")
print(f"\nAssuming perfect play, our team will earn {total_score} points.")
