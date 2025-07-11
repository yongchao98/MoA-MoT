# 1. Define the points for the winning team (X)
# Card points for the winning team in the optimal scenario.
# The last player (from the losing team) gives their tricks to the first player (from the winning team),
# allowing the winning team to capture all 100 points.
winning_team_card_points = 100

# Bet points for the winning team.
# One player on the winning team calls and succeeds in a Grand Tichu (+200).
winning_team_bet_points = 200

# Calculate the total score for the winning team, X.
X = winning_team_card_points + winning_team_bet_points

# 2. Define the points for the losing team (Y)
# Card points for the losing team in the optimal scenario.
losing_team_card_points = 0

# Bet points for the losing team.
# Both players on the losing team call Grand Tichu and fail (-200 each).
losing_team_bet_points = -200 + (-200)

# Calculate the total score for the losing team, Y.
Y = losing_team_card_points + losing_team_bet_points

# 3. Calculate the maximal difference X - Y
max_difference = X - Y

# 4. Print the final equation and the result.
print("This scenario maximizes the score difference (X - Y) under the given conditions.")
print("\nCalculating the winning team's score (X):")
print(f"Card Points = {winning_team_card_points}")
print(f"Bet Points = {winning_team_bet_points} (Successful Grand Tichu)")
print(f"X = {winning_team_card_points} + {winning_team_bet_points} = {X}")

print("\nCalculating the losing team's score (Y):")
print(f"Card Points = {losing_team_card_points}")
print(f"Bet Points = -200 - 200 = {losing_team_bet_points} (Two failed Grand Tichus)")
print(f"Y = {losing_team_card_points} + ({losing_team_bet_points}) = {Y}")

print("\nCalculating the maximal difference (X - Y):")
print(f"X - Y = {X} - ({Y}) = {max_difference}")