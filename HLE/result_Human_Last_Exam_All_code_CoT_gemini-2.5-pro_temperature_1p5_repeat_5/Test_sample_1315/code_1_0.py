import sys

# This script calculates the maximal score difference (X - Y) in a Tichu round
# under the specific condition that the winning team does not finish 1st and 2nd.

# Let Team A be the winning team (Score X) and Team B be the losing team (Score Y).
# The total score for a team is the sum of their card points and Tichu call points.
# We want to maximize X - Y.
# X - Y = (Card_Points_A + Tichu_Points_A) - (Card_Points_B + Tichu_Points_B)
# This can be rewritten as:
# X - Y = (Card_Points_A - Card_Points_B) + (Tichu_Points_A - Tichu_Points_B)

# --- Part 1: Maximizing the Card Point Difference ---
# To maximize the card point difference, Team A must get all positive points
# and Team B must get all negative points.
# Positive points: 4 Kings (40) + 4 Tens (40) + 4 Fives (20) + Dragon (25) = 125
# Negative points: Phoenix (-25)
max_card_points_A = 125
min_card_points_B = -25

# --- Part 2: Maximizing the Tichu Call Point Difference ---
# To maximize this difference, Team A succeeds in the highest call (Grand Tichu)
# and Team B fails the highest call.
# This happens if a player from Team A calls Grand Tichu and goes out first.
# Any other player who called automatically fails.
max_tichu_points_A = 200  # Successful Grand Tichu
min_tichu_points_B = -200 # Failed Grand Tichu

# This scenario requires a player from Team A to go out first, which satisfies the
# condition for the Grand Tichu bonus. The problem states the winning team (A)
# cannot finish 1st AND 2nd, which is not violated if a player from Team B
# finishes 2nd.

# --- Part 3: Calculate the final scores and the difference ---

# Calculate X, the maximal score for the winning team
X = max_card_points_A + max_tichu_points_A

# Calculate Y, the minimal score for the losing team
Y = min_card_points_B + min_tichu_points_B

# The maximal difference is X - Y
max_difference = X - Y

# --- Output the results step-by-step ---
print("To find the maximal value of X-Y, we must maximize X and minimize Y.")
print("-" * 60)

print("1. Calculating the maximal score for the winning team (X):")
print(f"   Maximal card points for winning team = {max_card_points_A}")
print(f"   Maximal Tichu bonus for winning team = {max_tichu_points_A}")
print(f"   X = {max_card_points_A} + {max_tichu_points_A} = {X}")
print()

print("2. Calculating the minimal score for the losing team (Y):")
print(f"   Minimal card points for losing team = {min_card_points_B}")
print(f"   Minimal Tichu penalty for losing team = {min_tichu_points_B}")
print(f"   Y = {min_card_points_B} + ({min_tichu_points_B}) = {Y}")
print()

print("3. Calculating the maximal difference (X - Y):")
print(f"   X - Y = {X} - ({Y}) = {max_difference}")

# Flushing the output buffer to ensure order.
sys.stdout.flush()
