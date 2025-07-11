# Plan: Calculate the maximal score difference (X-Y) in a Tichu round
# under the given constraints.

# 1. Maximize the winning team's (Team W) card points.
# The total points in the deck are 100. Team W can potentially capture all of them.
# This happens if the losing team's players either capture no points or if the
# last player to go out is from the losing team, giving their points to Team W.
winning_team_card_points = 100
losing_team_card_points = 100 - winning_team_card_points

# 2. Maximize the winning team's Tichu points.
# This is achieved by a player on Team W successfully calling and completing a Grand Tichu.
# A successful Grand Tichu grants 200 points.
winning_team_tichu_points = 200

# 3. Minimize the losing team's (Team L) Tichu points.
# This occurs when both players on Team L make calls and fail. The largest penalty
# is from one player failing a Grand Tichu (-200) and the other failing a Tichu (-100).
losing_team_failed_grand_tichu = -200
losing_team_failed_tichu = -100
losing_team_tichu_points = losing_team_failed_grand_tichu + losing_team_failed_tichu

# 4. Calculate the total scores for each team.
# X is the score of the winning team.
X = winning_team_card_points + winning_team_tichu_points

# Y is the score of the losing team.
Y = losing_team_card_points + losing_team_tichu_points

# 5. Calculate and print the final difference, X - Y.
# We explicitly show each number in the final equation.
difference = X - Y

print("This scenario assumes the following outcome for one round of Tichu:")
print("- A player from the winning team successfully calls 'Grand Tichu' (+200 points).")
print("- Both players from the losing team fail their 'Grand Tichu' and 'Tichu' calls (-200 and -100 points).")
print("- The winning team captures all 100 points from cards.")
print("-" * 30)
print(f"Winning Team Score (X):")
print(f"Card Points: {winning_team_card_points}")
print(f"Tichu Points: {winning_team_tichu_points}")
print(f"Total X = {winning_team_card_points} + {winning_team_tichu_points} = {X}")
print("-" * 30)
print(f"Losing Team Score (Y):")
print(f"Card Points: {losing_team_card_points}")
print(f"Tichu Points: {losing_team_failed_grand_tichu} + ({losing_team_failed_tichu}) = {losing_team_tichu_points}")
print(f"Total Y = {losing_team_card_points} + {losing_team_tichu_points} = {Y}")
print("-" * 30)
print("Maximal possible value of X - Y:")
print(f"{X} - ({Y}) = {difference}")
