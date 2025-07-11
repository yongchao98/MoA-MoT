#
# Plan to calculate the maximal possible value of X - Y in a Tichu round
# where X is the winning team's score and Y is the losing team's score,
# and the winning team does not go out first and second.
#

# Step 1: Maximize the score difference from card points.
# The total card points in a Tichu deck is 100.
# To maximize the difference, the winning team (Team A) gets all the points,
# and the losing team (Team B) gets none.
# This is possible if a player from Team B finishes last, giving all their
# point cards (from tricks and hand) to Team A.
winning_team_card_points = 100
losing_team_card_points = 0

# Step 2: Maximize the score difference from Tichu calls.
# Grand Tichu calls provide the largest point swings (+200 for success, -200 for failure).
#
# To maximize the winning team's bonus, one player makes a successful Grand Tichu call
# by going out first.
winning_team_call_bonus = 200  # One successful Grand Tichu

# To maximize the losing team's penalty, both players on the losing team
# attempt a Grand Tichu but fail (because a player from the winning team went out first).
losing_team_call_bonus = -200 + -200 # Two failed Grand Tichus

# Step 3: Calculate the total scores for the winning team (X) and the losing team (Y).
# The problem constraint (winning team does not finish 1-2) is satisfied
# with a finish order like: Winner1, Loser1, Winner2, Loser2.
X = winning_team_card_points + winning_team_call_bonus
Y = losing_team_card_points + losing_team_call_bonus

# Step 4: Calculate the final difference, X - Y.
difference = X - Y

# Print the results in a clear, step-by-step format.
print("To find the maximal possible value of X - Y, we construct the optimal scenario:")
print("\n1. Card Points Score:")
print(f"   - The winning team (Team A) captures all 100 card points.")
print(f"   - The losing team (Team B) captures 0 card points.")
print(f"   - Winning team card score = {winning_team_card_points}")
print(f"   - Losing team card score = {losing_team_card_points}")

print("\n2. Tichu Call Bonuses:")
print(f"   - A player on Team A makes a successful Grand Tichu call: +{winning_team_call_bonus} points.")
print(f"   - Both players on Team B make failed Grand Tichu calls: 2 * (-200) = {losing_team_call_bonus} points.")

print("\n3. Final Scores (X and Y):")
print(f"   - Winning Team's Score (X) = {winning_team_card_points} (cards) + {winning_team_call_bonus} (calls) = {X}")
print(f"   - Losing Team's Score (Y) = {losing_team_card_points} (cards) + ({losing_team_call_bonus}) (calls) = {Y}")

print("\n4. Maximal Difference (X - Y):")
print(f"   - The maximal value of X - Y = {X} - ({Y}) = {difference}")
