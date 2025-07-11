# Define scoring values for Tichu
grand_tichu_success = 200
grand_tichu_failure = -200
total_card_points_in_deck = 100

# Step 1: Maximize the card point difference.
# In the optimal scenario, the winning team (Team X) gets all card points.
# This happens if the last player to finish (from Team Y) collects all 100 points
# in their tricks and then must hand them to the first player out (from Team X).
winning_team_card_points = total_card_points_in_deck
losing_team_card_points = 0

# Step 2: Maximize the call point difference.
# To maximize the score difference, a player on the winning team makes a successful
# Grand Tichu call, while both players on the losing team fail their Grand Tichu calls.
winning_team_call_score = grand_tichu_success
losing_team_call_score = grand_tichu_failure + grand_tichu_failure # Two failed calls

# Step 3: Calculate the total scores for the winning team (X) and losing team (Y).
X = winning_team_card_points + winning_team_call_score
Y = losing_team_card_points + losing_team_call_score

# Step 4: Calculate the maximal difference, X - Y.
max_difference = X - Y

# Print the breakdown of the calculation as an equation.
print("Maximal Score Difference Calculation:")
print(f"Winning Team Score (X) = Card Points + Call Points")
print(f"X = {winning_team_card_points} + {winning_team_call_score} = {X}")
print("-" * 30)
print(f"Losing Team Score (Y) = Card Points + Call Points")
print(f"Y = {losing_team_card_points} + ({losing_team_call_score}) = {Y}")
print("-" * 30)
print("Final Difference (X - Y):")
# The final print statement shows each number in the final equation as requested.
print(f"{X} - ({Y}) = {max_difference}")
