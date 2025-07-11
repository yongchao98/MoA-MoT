# This script calculates the maximal score difference X-Y in a specific Tichu scenario.

# Step 1: Define the points for the winning team (Team A).
# A player from Team A calls a successful "Grand Tichu".
team_A_tichu_points = 200

# Team A manages to collect all 100 points from the cards.
team_A_card_points = 100

# The total score for Team A is X.
X = team_A_tichu_points + team_A_card_points

# Step 2: Define the points for the losing team (Team B).
# Both players from Team B call "Grand Tichu" and fail.
team_B_tichu_points = -200 + -200

# Team B collects no points from cards.
team_B_card_points = 0

# The total score for Team B is Y.
Y = team_B_tichu_points + team_B_card_points

# Step 3: Calculate the difference X - Y.
difference = X - Y

# Step 4: Print the detailed calculation and the final result.
print("This scenario calculates the maximal value of X-Y.")
print("The winning team's score, X, is composed of:")
print(f"  Tichu Points: A successful Grand Tichu = {team_A_tichu_points}")
print(f"  Card Points: All points captured = {team_A_card_points}")
print(f"X = {team_A_tichu_points} + {team_A_card_points} = {X}\n")

print("The losing team's score, Y, is composed of:")
print(f"  Tichu Points: Two failed Grand Tichus = {team_B_tichu_points}")
print(f"  Card Points: No points captured = {team_B_card_points}")
print(f"Y = {team_B_tichu_points} + {team_B_card_points} = {Y}\n")

print("The final difference is:")
print(f"X - Y = {X} - ({Y}) = {difference}")