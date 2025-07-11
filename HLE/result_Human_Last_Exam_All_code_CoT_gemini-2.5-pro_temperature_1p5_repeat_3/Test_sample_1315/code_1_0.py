# Tichu score components for the maximal difference scenario.

# For the winning team (Team A), we assume:
# - One successful "Grand Tichu" call.
# - Capturing all positive point cards (5s, 10s, Kings, Dragon).
team_A_tichu_score = 200
team_A_card_score = 125
X = team_A_tichu_score + team_A_card_score

# For the losing team (Team B), we assume:
# - One failed "Grand Tichu" call.
# - One failed "Tichu" call.
# - Capturing the trick with the Phoenix.
team_B_tichu_score = -200 + (-100)
team_B_card_score = -25
Y = team_B_tichu_score + team_B_card_score

# The maximal difference is X - Y.
max_difference = X - Y

# Print the final equation step by step.
print("To find the maximal value of X - Y, we construct the highest possible score for the winning team (X) and the lowest possible score for the losing team (Y).")
print("\n--- Calculating Winning Team's Score (X) ---")
print(f"Successful Grand Tichu: {team_A_tichu_score}")
print(f"Maximal Card Score (all positive cards): {team_A_card_score}")
print(f"X = {team_A_tichu_score} + {team_A_card_score} = {X}")


print("\n--- Calculating Losing Team's Score (Y) ---")
print(f"Failed Grand Tichu (-200) + Failed Tichu (-100): {team_B_tichu_score}")
print(f"Minimal Card Score (capturing the Phoenix): {team_B_card_score}")
print(f"Y = {team_B_tichu_score} + {team_B_card_score} = {Y}")

print("\n--- Calculating Maximal Difference (X - Y) ---")
print(f"X - Y = {X} - ({Y})")
print(f"X - Y = {max_difference}")