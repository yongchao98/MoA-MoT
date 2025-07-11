# Step 1: Calculate my meld points
run_points = 15
aces_100_points = 100
my_meld_total = run_points + aces_100_points

# Step 2: Calculate my partner's meld points
pinochle_points = 4
nines_of_trump_points = 2  # 1 point for each of the two 9s of trump
partner_meld_total = pinochle_points + nines_of_trump_points

# Step 3: Calculate the team's total meld points
team_meld_total = my_meld_total + partner_meld_total

# Step 4: Calculate the trick points
# In a 48-card pinochle deck, there are 8 Aces, 8 Tens, and 8 Kings.
# Each of these "counters" is worth 1 point when won in a trick.
counter_points = 8 + 8 + 8
# An additional point is awarded for winning the last trick.
last_trick_bonus = 1
# With all 8 aces and superior trump, our team can win every trick.
# Therefore, we win all possible trick points.
team_trick_total = counter_points + last_trick_bonus

# Step 5: Calculate the grand total for the hand
grand_total = team_meld_total + team_trick_total

# Final Output
print("Calculating the total points for the hand:")
print("\n--- Meld Points ---")
print(f"My Meld: {run_points} (Run) + {aces_100_points} (100 Aces) = {my_meld_total} points")
print(f"Partner's Meld: {pinochle_points} (Pinochle) + {nines_of_trump_points} (Two 9s of Trump) = {partner_meld_total} points")
print(f"Team Meld Total: {my_meld_total} + {partner_meld_total} = {team_meld_total} points")

print("\n--- Trick Points ---")
print("With all 8 aces and trump control, our team can win every trick.")
print(f"Total Counters (Aces, 10s, Kings): {counter_points} points")
print(f"Bonus for the last trick: {last_trick_bonus} point")
print(f"Team Trick Total: {counter_points} + {last_trick_bonus} = {team_trick_total} points")

print("\n--- Grand Total ---")
print(f"Final Score = Meld Points + Trick Points")
print(f"Final Score: {team_meld_total} + {team_trick_total} = {grand_total} points")

<<<146>>>