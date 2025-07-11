import sys

# This script calculates the maximal score difference (X-Y) in a Tichu round
# under the constraint that the winning team does not finish 1st and 2nd.

# Define the standard point values in Tichu
SUCCESSFUL_GRAND_TICHU = 200
FAILED_GRAND_TICHU = -200
SUCCESSFUL_TICHU = 100
FAILED_TICHU = -100
TOTAL_CARD_POINTS = 100

# To maximize X-Y, we must maximize X (winning team score) and minimize Y (losing team score).

# --- Step 1: Maximize the winning team's score (X) ---
# For Tichu calls, the best outcome is a successful Grand Tichu.
winning_team_tichu_points = SUCCESSFUL_GRAND_TICHU

# For card points, the best outcome is capturing all available points.
winning_team_card_points = TOTAL_CARD_POINTS

# --- Step 2: Minimize the losing team's score (Y) ---
# For Tichu calls, the worst outcome is for both players to fail their calls.
# To get the largest penalty, one player fails a Grand Tichu and the other fails a regular Tichu.
losing_team_failed_grand = FAILED_GRAND_TICHU
losing_team_failed_tichu = FAILED_TICHU
losing_team_tichu_points = losing_team_failed_grand + losing_team_failed_tichu

# For card points, the worst outcome is capturing no points.
losing_team_card_points = 0

# --- Step 3: Calculate the final scores X and Y ---
# A finishing order of Win(1st), Lose(2nd), Win(3rd), Lose(4th) meets the problem's constraint.
# In this scenario, all Tichu call outcomes and card point distributions are plausible.

# Calculate X (Winning Team Score)
X = winning_team_tichu_points + winning_team_card_points

# Calculate Y (Losing Team Score)
Y = losing_team_tichu_points + losing_team_card_points

# --- Step 4: Present the final calculation ---
# The prompt requires outputting each number in the final equation.

print("--- Maximizing the Score Difference (X-Y) ---")
print("This scenario assumes the following outcome:")
print("- Winning Team (A): One player calls and succeeds at Grand Tichu. The team collects all 100 card points.")
print("- Losing Team (B): One player calls and fails Grand Tichu, the other calls and fails a regular Tichu.")
print("\n--- Calculating the Scores ---")

print("\nWinning Team Score (X):")
print(f"  Tichu Points = {winning_team_tichu_points}")
print(f"  Card Points  = {winning_team_card_points}")
print(f"  X = {winning_team_tichu_points} + {winning_team_card_points} = {X}")

print("\nLosing Team Score (Y):")
print(f"  Tichu Points = {losing_team_failed_grand} + {losing_team_failed_tichu} = {losing_team_tichu_points}")
print(f"  Card Points  = {losing_team_card_points}")
print(f"  Y = {losing_team_tichu_points} + {losing_team_card_points} = {Y}")

print("\n--- Final Difference Calculation ---")
print("The maximal possible value of X-Y is calculated as:")
print(f"({X}) - ({Y}) = {X - Y}")
